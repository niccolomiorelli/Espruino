/*
 * This file is part of Espruino, a JavaScript interpreter for Microcontrollers
 *
 * Copyright (C) 2021 Gordon Williams <gw@pur3.co.uk>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * ----------------------------------------------------------------------------
 * Heart rate monitoring — TRUST-PPG4 algorithm
 *
 * Replaces heartrate_vc31_binary.c for Bangle.js 2 (NRF52840).
 * Implements the same external API:
 *   hrm_init(), hrm_new(), hrm_get_hrm_info(), hrm_get_hrm_raw_info()
 *
 * Pipeline:
 *  1. Input: PPG sample + ACC sample (last known ACC used when ACC rate < PPG rate).
 *  2. TRUST-PPG4: bandpass + z-score + NLMS adaptive filter + FFT spectral
 *     peak detection + dual-measure Kalman filter for motion-robust HR.
 *  3. Output: hrmInfo.bpm10 (bpm×10), hrmInfo.confidence (0–100).
 * ----------------------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include "heartrate.h"
#include "hrm.h"
#include "jshardware.h"
#include "trust_ppg4/trust_ppg4_heartrate.h"

/* =========================================================
 * Global state required by the Espruino HRM subsystem
 * ========================================================= */

HrmInfo    hrmInfo;
HrmSample  hrmSamples[HRMSAMPLE_MAX];
uint8_t    hrmSampleCount;

/* Last known ACC values (used when ACC rate < PPG rate) */
static int16_t last_ax, last_ay, last_az;

/* Last reported BPM to avoid duplicate events */
static int  lastHR_bpm10;

/* ---- hrm_init ---- */
void hrm_init(void)
{
    memset(&hrmInfo, 0, sizeof(hrmInfo));
    /* lastBeatTime tracks the time of the last call to hrm_new() */
    hrmInfo.lastBeatTime = jshGetSystemTime();

    trust_ppg4_heartrate_init();

    last_ax = last_ay = last_az = 0;
    lastHR_bpm10   = 0;
}

/* ---- hrm_new ---- */
bool hrm_new(int hrmValue, Vector3 *acc)
{
    /* Time since last sample */
    JsSysTime time = jshGetSystemTime();
    int timeDiff = (int)(jshGetMillisecondsFromTime(time - hrmInfo.lastBeatTime) + 0.5);
    hrmInfo.lastBeatTime = time;
    if (timeDiff < 1)    timeDiff = 1;
    if (timeDiff > 2000) timeDiff = 2000;

    /* Update raw/avg/filtered for HRM-raw events */
    if (hrmValue < HRMVALUE_MIN) hrmValue = HRMVALUE_MIN;
    if (hrmValue > HRMVALUE_MAX) hrmValue = HRMVALUE_MAX;
    hrmInfo.raw = (HrmValueType)hrmValue;
    if (hrmPollInterval > 30) /* 40 ms = 25 Hz */
        hrmInfo.avg = (int16_t)(((int)hrmInfo.avg * 7  + hrmValue) >> 3);
    else                       /* 20 ms = 50 Hz */
        hrmInfo.avg = (int16_t)(((int)hrmInfo.avg * 15 + hrmValue) >> 4);
    hrmInfo.filtered = (int16_t)(hrmValue - hrmInfo.avg);

    /* Update ACC — keep last known value if acc is NULL */
    if (acc) {
        last_ax = acc->x;
        last_ay = acc->y;
        last_az = acc->z;
    }

    int16_t dt = (int16_t)(timeDiff < 32767 ? timeDiff : 32767);
    int new_bpm10 = trust_ppg4_heartrate(dt,
                                          (ppg_t)hrmValue,
                                          (accel_t)last_ax,
                                          (accel_t)last_ay,
                                          (accel_t)last_az);
    int new_conf = trust_ppg4_get_last_confidence(); /* 0–100 */

    hrmInfo.isBeat = false;
    bool hadBeat = false;

    if (new_bpm10 > 0 && new_bpm10 != lastHR_bpm10) {
        hrmInfo.bpm10      = (uint16_t)new_bpm10;
        hrmInfo.confidence = (uint8_t)(new_conf > 100 ? 100 : new_conf);
        lastHR_bpm10       = new_bpm10;
        hrmInfo.isBeat     = true;
        hadBeat            = true;
    }

    /* Fill sample buffer for HRM-raw events */
    if (hrmSampleCount < HRMSAMPLE_MAX) {
        HrmSample *sample = &hrmSamples[hrmSampleCount++];
        sample->bpm10      = hrmInfo.bpm10;
        sample->confidence = hrmInfo.confidence;
        sample->raw        = hrmInfo.raw;
        sample->avg        = hrmInfo.avg;
        sample->filtered   = hrmInfo.filtered;
    }

    return hadBeat;
}

/* ---- hrm_get_hrm_info ---- */
void hrm_get_hrm_info(JsVar *o)
{
    /* TRUST-PPG4 is window-based, no individual beat timestamps.
       Emit an empty history array to keep the JS interface consistent. */
    JsVar *a = jsvNewEmptyArray();
    if (a) {
        jsvObjectSetChildAndUnLock(o, "history", a);
    }
}

/* ---- hrm_get_hrm_raw_info ---- */
void hrm_get_hrm_raw_info(JsVar *o)
{
    jsvObjectSetBoolChild(o, "isBeat", hrmInfo.isBeat);
    uint32_t sc, wc;
    trust_ppg4_get_timing(&sc, &wc);
    jsvObjectSetIntegerChild(o, "cyclesSample", (JsVarInt)sc);
    jsvObjectSetIntegerChild(o, "cyclesWindow",  (JsVarInt)wc);
}

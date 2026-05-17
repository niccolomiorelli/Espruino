#ifndef TRUST_PPG4_HEARTRATE_H
#define TRUST_PPG4_HEARTRATE_H

#include "trust_ppg4_types.h"

/** Reset all internal state. Must be called before the first hrm sample. */
void trust_ppg4_heartrate_init(void);

/**
 * Process one PPG+ACC sample (with linear interpolation for missing samples).
 * @param delta_ms   Time in ms since the previous call.
 * @param ppg        Raw PPG sensor value.
 * @param accx/y/z   Raw accelerometer values (any consistent unit).
 * @return           Heart rate estimate in bpm*10 (e.g. 720 = 72.0 bpm),
 *                   or 0 if no estimate is available yet.
 */
int trust_ppg4_heartrate(time_delta_ms_t delta_ms,
                          ppg_t ppg,
                          accel_t accx, accel_t accy, accel_t accz);

/**
 * Return the confidence of the last spectral window as a percentage (0-100).
 * Valid only after the first window has been computed.
 */
int trust_ppg4_get_last_confidence(void);

#endif

#ifndef KALMAN_FILTER_FLOAT_H_
#define KALMAN_FILTER_FLOAT_H_

#include "trust_ppg4_types.h"

/* Motion-state-based process noise */
float find_Qk_float2(int state);
/* Confidence-based measurement noise */
float find_Rk_float(float conf);
float find_Rk_float2(float conf);

float clampf(float x, float lo, float hi);

void kalman_HR_init_float(float *HR_est, float *P);
void kalman_HR_2meas_init_float(float *HR_est, float *P);

/* Single-measurement Kalman update */
void kalman_HR_estimation_float2(float HR_meas, int state, float conf,
                                  float *HR_est, float *P);

/* Dual-measurement (PPG + ACC) Kalman update */
void kalman_HR_estimation_2measures_float2(float HR_meas_PPG, float HR_meas_ACC,
                                            int state, float conf,
                                            float *HR_est, float *P);

#endif

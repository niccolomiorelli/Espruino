#include "kalman_filter_float.h"
#include <math.h>

#define KF_DELTA_MIN   0.05f
#define KF_DELTA_MAX   0.83f
#define KF_Z_SCORE     1.96f
#define KF_MIN_CONF    0.001f

#define STATE_STATIONARY    0
#define STATE_TRANSITION    1
#define STATE_MOTION_LOW    2
#define STATE_MOTION_MEDIUM 3
#define STATE_MOTION_HIGH   4

float clampf(float x, float lo, float hi)
{
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

/* Process noise: function of motion state */
float find_Qk_float2(int state)
{
    float var;
    if (state == STATE_STATIONARY || state == STATE_TRANSITION || state == STATE_MOTION_LOW)
        var = 0.01f;
    else if (state == STATE_MOTION_MEDIUM)
        var = 0.02f;
    else
        var = 0.03f;
    return var * var;
}

/* Measurement noise: function of confidence (wide range) */
float find_Rk_float(float conf)
{
    if (conf < KF_MIN_CONF) conf = KF_MIN_CONF;
    float delta_Hz = KF_DELTA_MIN + (1.0f - conf) * (KF_DELTA_MAX - KF_DELTA_MIN);
    float sigma = delta_Hz / KF_Z_SCORE;
    return sigma * sigma;
}

/* Measurement noise: function of confidence (tighter range) */
float find_Rk_float2(float conf)
{
    if (conf < KF_MIN_CONF) conf = KF_MIN_CONF;
    float delta_Hz = KF_DELTA_MIN + (1.0f - conf) * (1.2f - KF_DELTA_MIN);
    float sigma = delta_Hz / KF_Z_SCORE;
    return sigma * sigma;
}

void kalman_HR_init_float(float *HR_est, float *P)
{
    *HR_est = 1.0f;
    *P      = 0.2f;
}

void kalman_HR_2meas_init_float(float *HR_est, float *P)
{
    *HR_est = 1.0f;
    *P      = 0.1f;
}

/* Single-measurement Kalman update */
void kalman_HR_estimation_float2(float HR_meas, int state, float conf,
                                  float *HR_est, float *P)
{
    float Qk = find_Qk_float2(state);
    float Rk = find_Rk_float(conf);

    float HR_pred = *HR_est;
    float P_pred  = *P + Qk;

    if (conf > 0.3f) {
        float K  = P_pred / (P_pred + Rk);
        *HR_est  = HR_pred + K * (HR_meas - HR_pred);
        *P       = (1.0f - K) * P_pred;
    } else {
        *HR_est = HR_pred;
        *P      = P_pred;
    }
}

/* Dual-measurement (PPG + ACC-model) Kalman update */
void kalman_HR_estimation_2measures_float2(float HR_meas_PPG, float HR_meas_ACC,
                                            int state, float conf,
                                            float *HR_est, float *P)
{
    float Qk     = find_Qk_float2(state);
    float Rk_PPG = find_Rk_float(conf);
    float Rk_ACC = (0.25f / KF_Z_SCORE) * (0.25f / KF_Z_SCORE);

    float HR_pred = *HR_est;
    float P_pred  = *P + Qk;

    if (conf >= 0.3f) {
        float det = P_pred * (Rk_PPG + Rk_ACC) + Rk_PPG * Rk_ACC;
        float K1  = P_pred * Rk_ACC / det;
        float K2  = P_pred * Rk_PPG / det;

        float y1 = HR_meas_PPG - HR_pred;
        float y2 = HR_meas_ACC - HR_pred;

        *HR_est = HR_pred + K1 * y1 + K2 * y2;
        *P      = (1.0f - (K1 + K2)) * P_pred;
    } else {
        *HR_est = HR_pred;
        *P      = P_pred;
    }
}

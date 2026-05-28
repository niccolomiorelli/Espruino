/* ----------------------------------------------------
 * TRUST-PPG 4 — adapted for Espruino / Bangle.js
 *
 * Pipeline:
 *  1. Bandpass filter + z-score standardisation of PPG and ACC.
 *  2. Multi-input NLMS adaptive filter for motion-artefact removal.
 *  3. State machine: Stationary / Transition / Low / Medium / High motion.
 *  4. FFT-based spectral analysis (PPG or NLMS output, depending on state).
 *  5. Peak selection with motion-state-dependent frequency range.
 *  6. Confidence: C1 (peak/total power), C2 (peak/local power), C3 (consistency).
 *  7. Dual-measure Kalman filter (PPG peak + ACC-model prior) for HR tracking.
 *  8. Output: HR in bpm*10.
 * ---------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "trust_ppg4_types.h"
#include "bandpass_filter.h"
#include "rolling_stats_float.h"
#include "fft_library_float.h"
#include "adaptive_filter_float.h"
#include "kalman_filter_float.h"
#include "trust_ppg4_heartrate.h"

/* ---- DWT cycle counter (nRF52840 only) ---- */
#ifdef NRF52840_XXAA
#  define TP4_DWT_CTRL    (*((volatile uint32_t*)0xE0001000u))
#  define TP4_DWT_CYCCNT  (*((volatile uint32_t*)0xE0001004u))
#  define TP4_DEMCR       (*((volatile uint32_t*)0xE000EDFCu))
#  define TP4_DWT_ENABLE() do { \
       TP4_DEMCR     |= (1u << 24); \
       TP4_DWT_CYCCNT = 0u;          \
       TP4_DWT_CTRL  |= 1u;          \
   } while(0)
#else
#  define TP4_DWT_ENABLE()  do {} while(0)
#  define TP4_DWT_CYCCNT    0u
#endif

/** Cycles consumed by one 25 Hz call (per-sample cost). */
static volatile uint32_t tp4_cycles_sample = 0;
/** Cycles consumed by the FFT window block (every 64 calls). */
static volatile uint32_t tp4_cycles_window = 0;

/* ---- Parameters ---- */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef TP4_WINDOW_LEN
#define TP4_WINDOW_LEN   128
#endif
#define TP4_WINDOW_STEP   64
#define TP4_SAMPLING_FREQ 25
#ifndef TP4_N_PAD
#define TP4_N_PAD        512
#endif

/* Physiological HR range (36–210 bpm) as FFT bin indices @25 Hz / 512-point */
#define LOW_FREQ_PHY_I   13 //25
#define HIGH_FREQ_PHY_I  73 //143
#define DELTA_HR_RANGE   17 //34   /* ~50 bpm */

/* State machine */
#define THRESHOLD_MOTION1   200.0f
#define THRESHOLD_MOTION2   300.0f
#define THRESHOLD_MOTION3  2000.0f
#define STATE_STATIONARY    0
#define STATE_TRANSITION    1
#define STATE_MOTION_LOW    2
#define STATE_MOTION_MEDIUM 3
#define STATE_MOTION_HIGH   4
#define N_CONSEC_STATIONARY_TO_MOTION 4
#define N_CONSEC_OTHER_CHANGES        2

/* NLMS */
#define MU 0.2f

/* Confidence */
#define DELTA_PEAK_i    16
#define DELTA_C3_i       8
#define SCALE_C1        10.0f
#define SCALE_C2        10.0f
#define MIDPOINT_C1     0.35f
#define MIDPOINT_C2     0.65f
#define ALPHA_W         0.45f
#define BETA_W          0.45f
#define GAMMA_W         0.10f

#define N_LAST_SAVED    3
#define EPSILON         1e-10f

/* ---- Static state ---- */

static ppg_t   ppg_last;
static accel_t accx_last, accy_last, accz_last;

/* Bandpass filters */
static BPFilter bpFilter_ppg;
static BPFilter bpFilter_accx;
static BPFilter bpFilter_accy;
static BPFilter bpFilter_accz;

/* Rolling statistics */
static Stats_float stats_ppg;
static Stats_float stats_accx;
static Stats_float stats_accy;
static Stats_float stats_accz;
static Stats_float stats_acc;
static Stats_float stats_acc_raw;

/* NLMS adaptive filter */
static AdaptFilter_multi_float NLMS_filter;

/* Circular sample buffers */
static float   signal_PPG_buffer[TP4_WINDOW_LEN];
static float   signal_NLMS_buffer[TP4_WINDOW_LEN];
static accel_t acc_magnitude[TP4_WINDOW_LEN];
static int     signal_buffer_next_i;
static int     samples_since_last_HR;

/* FFT workspace */
static float              windowed_signal_out[TP4_WINDOW_LEN]; /* also used as input to Hann window (in-place) */
static complex_number_float fft_input_padded[TP4_N_PAD];        /* first WINDOW_LEN slots used before zero-padding */

/* State machine */
static int state;
static int pending_state;
static int pending_count;

/* ACC→HR linear model (b is updated online) */
static float model_a;
static float model_b;

/* History of last N spectral peaks and their confidences */
static int   last_peaks_i[N_LAST_SAVED];
static float last_peaks[N_LAST_SAVED];
static float last_confs[N_LAST_SAVED];

/* Kalman filters */
static float HR_freq_est2;
static float P_kf1;
static float HR_freq_est3;
static float P_kf2;

static bool  high_confidence;
static bool  first_window;

static int HR; /* current HR estimate in bpm*10 */

/* ---- Helper functions ---- */

static float HR_from_ACC(float rms, float a, float b)
{
    float out = a * rms + b;
    if (out > 210.0f) out = 210.0f;
    if (out <  36.0f) out =  36.0f;
    return out;
}

static void update_ACC_model(float HR_PPG, float HR_ACC, float *a, float *b, int st)
{
    float lr = 0.25f;
    (void)a; /* slope kept fixed */
    *b += lr * (HR_PPG - HR_ACC);
    if (*b < -15.0f) *b = -15.0f;
    if (*b >  70.0f) *b =  70.0f;
    if (st == STATE_MOTION_HIGH   && *b <  5.0f) *b =  5.0f;
}

static float sigmoid_conf(float x, float midpoint, float scale)
{
    return 1.0f / (1.0f + expf(-scale * (x - midpoint)));
}

static void weighted_avg_conf(const float *peaks, const float *confs, int len,
                               float *wavg, float *c_est)
{
    float sw  = 0.0f, swp = 0.0f, sw2 = 0.0f;
    for (int i = 0; i < len; i++) {
        swp += peaks[i] * confs[i];
        sw  += confs[i];
        sw2 += confs[i] * confs[i];
    }
    sw += EPSILON;
    *wavg  = swp / sw;
    *c_est = sw2 / sw;
}

/* ---- Init ---- */

void trust_ppg4_heartrate_init(void)
{
    ppg_last = accx_last = accy_last = accz_last = 0;

    BPFilter_init(&bpFilter_ppg);
    BPFilter_init(&bpFilter_accx);
    BPFilter_init(&bpFilter_accy);
    BPFilter_init(&bpFilter_accz);

    rolling_stats_reset_float(&stats_ppg);
    rolling_stats_reset_float(&stats_accx);
    rolling_stats_reset_float(&stats_accy);
    rolling_stats_reset_float(&stats_accz);
    rolling_stats_reset_float(&stats_acc);
    rolling_stats_reset_float(&stats_acc_raw);

    AdaptFilter_multi_init_float(&NLMS_filter);

    memset(signal_PPG_buffer,  0, sizeof(signal_PPG_buffer));
    memset(signal_NLMS_buffer, 0, sizeof(signal_NLMS_buffer));
    memset(acc_magnitude,      0, sizeof(acc_magnitude));
    memset(windowed_signal_out, 0, sizeof(windowed_signal_out));

    signal_buffer_next_i  = 0;
    samples_since_last_HR = 0;
    first_window          = true;

    state         = STATE_STATIONARY;
    pending_state = -1;
    pending_count = 0;

    model_a = 0.0093f;
    model_b = 14.96f;

    for (int i = 0; i < N_LAST_SAVED; i++) {
        last_peaks_i[i] = 0;
        last_peaks[i]   = 0.0f;
        last_confs[i]   = 0.0f;
    }

    kalman_HR_init_float(&HR_freq_est2, &P_kf1);
    kalman_HR_2meas_init_float(&HR_freq_est3, &P_kf2);

    high_confidence = false;
    HR = 0;
}

/* ---- Confidence getter ---- */

int trust_ppg4_get_last_confidence(void)
{
    float c = last_confs[0];
    if (c < 0.0f) c = 0.0f;
    if (c > 1.0f) c = 1.0f;
    return (int)(c * 100.0f);
}

/* ---- Core algorithm (one 25 Hz sample) ---- */

static int main_algorithm_trust_ppg4(time_delta_ms_t delta_ms,
                                      ppg_t ppg,
                                      accel_t accx, accel_t accy, accel_t accz)
{
    (void)delta_ms;
    TP4_DWT_ENABLE();

    /* --- Bandpass filter --- */
    accel_t acc_raw = (accel_t)sqrtf((float)accx*(float)accx +
                                      (float)accy*(float)accy +
                                      (float)accz*(float)accz);

    BPFilter_put(&bpFilter_ppg, ppg);
    ppg_t ppg_f = (ppg_t)BPFilter_get(&bpFilter_ppg);

    BPFilter_put(&bpFilter_accx, (ppg_t)accx);
    accel_t accx_f = (accel_t)BPFilter_get(&bpFilter_accx);
    BPFilter_put(&bpFilter_accy, (ppg_t)accy);
    accel_t accy_f = (accel_t)BPFilter_get(&bpFilter_accy);
    BPFilter_put(&bpFilter_accz, (ppg_t)accz);
    accel_t accz_f = (accel_t)BPFilter_get(&bpFilter_accz);

    accel_t acc_f = (accel_t)sqrtf((float)accx_f*(float)accx_f +
                                    (float)accy_f*(float)accy_f +
                                    (float)accz_f*(float)accz_f);

    /* --- Z-score standardisation --- */
    rolling_stats_addValue_float((float)ppg_f, &stats_ppg);
    if (rolling_stats_get_variance_float(&stats_ppg) < 0.0f)
        rolling_stats_reset_float(&stats_ppg);
    float std_ppg = rolling_stats_get_standard_deviation_float(&stats_ppg);
    if (std_ppg == 0.0f) std_ppg = 1.0f;
    float ppg_std = ((float)ppg_f - rolling_stats_get_mean_float(&stats_ppg)) / std_ppg;

    rolling_stats_addValue_float((float)acc_f, &stats_acc);
    if (rolling_stats_get_variance_float(&stats_acc) < 0.0f)
        rolling_stats_reset_float(&stats_acc);
    float std_acc = rolling_stats_get_standard_deviation_float(&stats_acc);
    if (std_acc == 0.0f) std_acc = 1.0f;
    float acc_std = ((float)acc_f - rolling_stats_get_mean_float(&stats_acc)) / std_acc;
    (void)acc_std; /* used indirectly below */

    rolling_stats_addValue_float((float)accx_f, &stats_accx);
    if (rolling_stats_get_variance_float(&stats_accx) < 0.0f)
        rolling_stats_reset_float(&stats_accx);
    float accx_std = ((float)accx_f - rolling_stats_get_mean_float(&stats_accx)) / std_acc;

    rolling_stats_addValue_float((float)accy_f, &stats_accy);
    if (rolling_stats_get_variance_float(&stats_accy) < 0.0f)
        rolling_stats_reset_float(&stats_accy);
    float accy_std = ((float)accy_f - rolling_stats_get_mean_float(&stats_accy)) / std_acc;

    rolling_stats_addValue_float((float)accz_f, &stats_accz);
    if (rolling_stats_get_variance_float(&stats_accz) < 0.0f)
        rolling_stats_reset_float(&stats_accz);
    float accz_std = ((float)accz_f - rolling_stats_get_mean_float(&stats_accz)) / std_acc;

    rolling_stats_addValue_float((float)acc_raw, &stats_acc_raw);
    if (rolling_stats_get_variance_float(&stats_acc_raw) < 0.0f)
        rolling_stats_reset_float(&stats_acc_raw);
    float std_acc_raw = rolling_stats_get_standard_deviation_float(&stats_acc_raw);
    if (std_acc_raw == 0.0f) std_acc_raw = 1.0f;

    /* --- NLMS adaptive filter for motion-artefact removal --- */
    AdaptFilter_multi_put_float(&NLMS_filter, accx_std, accy_std, accz_std, ppg_std);
    float out_NLMS = AdaptFilter_multi_get_float(&NLMS_filter, MU);

    /* --- Fill circular buffers --- */
    signal_PPG_buffer[signal_buffer_next_i]  = ppg_std;
    signal_NLMS_buffer[signal_buffer_next_i] = out_NLMS;
    acc_magnitude[signal_buffer_next_i]      = acc_raw;
    signal_buffer_next_i = (signal_buffer_next_i + 1) % TP4_WINDOW_LEN;

    samples_since_last_HR++;

    if (samples_since_last_HR == TP4_WINDOW_STEP) {
        if (first_window) {
            first_window = false;
            samples_since_last_HR = 0;
            tp4_cycles_sample = TP4_DWT_CYCCNT;
            return HR;
        }

        samples_since_last_HR = 0;
        uint32_t win_start = TP4_DWT_CYCCNT;

        /* --- State machine --- */
        int desired;
        if      (std_acc_raw < THRESHOLD_MOTION1) desired = STATE_STATIONARY;
        else if (std_acc_raw < THRESHOLD_MOTION2) desired = STATE_MOTION_LOW;
        else if (std_acc_raw < THRESHOLD_MOTION3) desired = STATE_MOTION_MEDIUM;
        else                                      desired = STATE_MOTION_HIGH;

        if (state == STATE_STATIONARY) {
            if (desired != STATE_STATIONARY) {
                state         = STATE_TRANSITION;
                pending_state = desired;
                pending_count = 1;
            } else {
                pending_state = -1;
                pending_count = 0;
            }
        } else if (state == STATE_TRANSITION) {
            if (desired == STATE_STATIONARY) {
                state         = STATE_STATIONARY;
                pending_state = -1;
                pending_count = 0;
            } else {
                pending_count++;
                pending_state = desired;
                if (pending_count >= N_CONSEC_STATIONARY_TO_MOTION) {
                    state         = pending_state;
                    pending_state = -1;
                    pending_count = 0;
                }
            }
        } else {
            if (desired == state) {
                pending_state = -1;
                pending_count = 0;
            } else {
                if (pending_state != desired) {
                    pending_state = desired;
                    pending_count = 1;
                } else {
                    pending_count++;
                    if (pending_count >= N_CONSEC_OTHER_CHANGES) {
                        state         = pending_state;
                        pending_state = -1;
                        pending_count = 0;
                    }
                }
            }
        }

        /* --- RMS of ACC magnitude over window --- */
        float cumul = 0.0f;
        for (int i = 0; i < TP4_WINDOW_LEN; i++) {
            float v = (float)acc_magnitude[i];
            cumul += v * v;
        }
        float rms_acc = sqrtf(cumul / TP4_WINDOW_LEN);

        /* --- ACC-model HR prior --- */
        float HR_ACC = HR_from_ACC(rms_acc, model_a, model_b);

        /* --- Frequency search range (FFT bin indices) --- */
        int lower_limit  = LOW_FREQ_PHY_I;
        int higher_limit = HIGH_FREQ_PHY_I;
        if (state == STATE_MOTION_MEDIUM) {
            lower_limit  = (int)ceilf( 60.0f * TP4_N_PAD / (TP4_SAMPLING_FREQ * 60.0f));
            higher_limit = (int)floorf(135.0f * TP4_N_PAD / (TP4_SAMPLING_FREQ * 60.0f));
        } else if (state == STATE_MOTION_HIGH) {
            lower_limit  = (int)ceilf( 100.0f * TP4_N_PAD / (TP4_SAMPLING_FREQ * 60.0f));
            higher_limit = (int)floorf(180.0f * TP4_N_PAD / (TP4_SAMPLING_FREQ * 60.0f));
        }

        /* --- FFT --- */
        const float *src = (state == STATE_STATIONARY) ? signal_PPG_buffer
                                                        : signal_NLMS_buffer;
        for (int i = 0; i < TP4_WINDOW_LEN; i++) {
            int bi = buffer_index_plus_fftLib(signal_buffer_next_i, i, TP4_WINDOW_LEN);
            windowed_signal_out[i] = src[bi];
        }
        apply_hann_window_fftLib_float(windowed_signal_out, windowed_signal_out, TP4_WINDOW_LEN);
        for (int i = 0; i < TP4_WINDOW_LEN; i++) {
            fft_input_padded[i].real = windowed_signal_out[i];
            fft_input_padded[i].imag = 0.0f;
        }
        zero_pad_fftLib_float(fft_input_padded, fft_input_padded, TP4_WINDOW_LEN);
        FFT_fftLib_float(fft_input_padded, 1.0f);

        /* --- Spectral magnitude + total power --- */
        float fft_magnitude[HIGH_FREQ_PHY_I + 1];
        float total_power = 0.0f;
        for (int i = 0; i <= HIGH_FREQ_PHY_I; i++) {
            if (i < LOW_FREQ_PHY_I) {
                fft_magnitude[i] = 0.0f;
            } else {
                fft_magnitude[i] = (fft_input_padded[i].real * fft_input_padded[i].real +
                                    fft_input_padded[i].imag * fft_input_padded[i].imag)
                                   / TP4_N_PAD;
                total_power += fft_magnitude[i];
            }
        }

        /* --- Peak detection --- */
        float max_mag = 0.0f;
        int   dom_idx = 0;
        for (int i = lower_limit + 1; i < higher_limit; i++) {
            if (fft_magnitude[i] > fft_magnitude[i - 1] &&
                fft_magnitude[i] > fft_magnitude[i + 1] &&
                fft_magnitude[i] > max_mag) {
                max_mag = fft_magnitude[i];
                dom_idx = i;
            }
        }

        bool peak_found = (dom_idx != 0);
        if (!peak_found) dom_idx = last_peaks_i[0];

        /* Shift history */
        for (int i = N_LAST_SAVED - 1; i > 0; i--) {
            last_peaks_i[i] = last_peaks_i[i - 1];
            last_peaks[i]   = last_peaks[i - 1];
        }
        last_peaks_i[0] = dom_idx;
        last_peaks[0]   = (float)dom_idx * TP4_SAMPLING_FREQ / (float)TP4_N_PAD;

        /* --- Confidence coefficients --- */
        /* C1: peak power / total power */
        float peak_power = 0.0f;
        int pl = dom_idx - DELTA_PEAK_i / 2;
        int ph = dom_idx + DELTA_PEAK_i / 2;
        if (pl < LOW_FREQ_PHY_I) pl = LOW_FREQ_PHY_I;
        if (ph > HIGH_FREQ_PHY_I) ph = HIGH_FREQ_PHY_I;
        for (int i = pl; i <= ph; i++) peak_power += fft_magnitude[i];
        float c1 = peak_power / (total_power + EPSILON);

        /* C2: peak power / local window power (~50 bpm) */
        float peak_power2 = peak_power; /* same window */
        float local_power = 0.0f;
        int hl = dom_idx - DELTA_HR_RANGE / 2;
        int hh = dom_idx + DELTA_HR_RANGE / 2;
        if (hl < LOW_FREQ_PHY_I) hl = LOW_FREQ_PHY_I;
        if (hh > HIGH_FREQ_PHY_I) hh = HIGH_FREQ_PHY_I;
        for (int i = hl; i <= hh; i++) local_power += fft_magnitude[i];
        float c2 = (local_power > 0.0f) ? (peak_power2 / local_power) : 0.0f;

        /* C3: temporal consistency */
        float c3 = 0.0f;
        if (abs(last_peaks_i[0] - last_peaks_i[1]) <= DELTA_C3_i / 2) {
            c3 = 0.75f;
            if (abs(last_peaks_i[0] - last_peaks_i[2]) <= DELTA_C3_i / 2)
                c3 = 1.0f;
        }

        float c1s = sigmoid_conf(c1, MIDPOINT_C1, SCALE_C1);
        float c2s = sigmoid_conf(c2, MIDPOINT_C2, SCALE_C2);
        float conf = ALPHA_W * c1s + BETA_W * c2s + GAMMA_W * c3;

        if (state == STATE_TRANSITION) conf -= 0.4f;
        if (!peak_found)               conf  = 0.0f;
        if (conf < 0.0f) conf = 0.0f;
        if (conf > 1.0f) conf = 1.0f;

        /* Shift confidence history */
        for (int i = N_LAST_SAVED - 1; i > 0; i--)
            last_confs[i] = last_confs[i - 1];
        last_confs[0] = conf;

        /* --- Kalman tracking --- */
        float HR_meas = last_peaks[0]; /* Hz */

        kalman_HR_estimation_float2(HR_meas, state, conf, &HR_freq_est2, &P_kf1);
        kalman_HR_estimation_2measures_float2(HR_meas, HR_ACC / 60.0f,
                                              state, conf,
                                              &HR_freq_est3, &P_kf2);

        HR = (int)(HR_freq_est3 * 600.0f); /* bpm*10 */

        /* --- High-confidence check --- */
        int cnt = 0;
        for (int i = 0; i < N_LAST_SAVED; i++)
            if (last_confs[i] >= 0.8f) cnt++;
        high_confidence = (cnt >= 2);

        /* --- Update ACC model when confidence is high --- */
        if (high_confidence) {
            float HR_avg, c_avg;
            weighted_avg_conf(last_peaks, last_confs, N_LAST_SAVED, &HR_avg, &c_avg);
            update_ACC_model(HR_avg * 60.0f, HR_ACC, &model_a, &model_b, state);
        }
        tp4_cycles_window = TP4_DWT_CYCCNT - win_start;
    }

    tp4_cycles_sample = TP4_DWT_CYCCNT;
    return HR;
}

/* ---- Timing getter ---- */
void trust_ppg4_get_timing(uint32_t *sample_cycles, uint32_t *window_cycles) {
    *sample_cycles = tp4_cycles_sample;
    *window_cycles = tp4_cycles_window;
}

/* ---- Public entry point with linear interpolation ---- */

int trust_ppg4_heartrate(time_delta_ms_t delta_ms,
                          ppg_t ppg,
                          accel_t accx, accel_t accy, accel_t accz)
{
    int final_HR;

    if (delta_ms < 60) {
        final_HR = main_algorithm_trust_ppg4(delta_ms, ppg, accx, accy, accz);
    } else if (delta_ms < 100) {
        int dm2   = delta_ms / 2;
        ppg_t pi  = (ppg_t)((ppg + ppg_last) / 2);
        accel_t xi = (accel_t)((accx + accx_last) / 2);
        accel_t yi = (accel_t)((accy + accy_last) / 2);
        accel_t zi = (accel_t)((accz + accz_last) / 2);
        main_algorithm_trust_ppg4(dm2, pi, xi, yi, zi);
        final_HR = main_algorithm_trust_ppg4(delta_ms, ppg, accx, accy, accz);
    } else {
        int dm1   = delta_ms / 3;
        int dm2   = delta_ms * 2 / 3;
        ppg_t pi1 = (ppg_t)(ppg / 3 + ppg_last * 2 / 3);
        ppg_t pi2 = (ppg_t)(ppg * 2 / 3 + ppg_last / 3);
        accel_t xi1 = (accel_t)(accx / 3 + accx_last * 2 / 3);
        accel_t yi1 = (accel_t)(accy / 3 + accy_last * 2 / 3);
        accel_t zi1 = (accel_t)(accz / 3 + accz_last * 2 / 3);
        accel_t xi2 = (accel_t)(accx * 2 / 3 + accx_last / 3);
        accel_t yi2 = (accel_t)(accy * 2 / 3 + accy_last / 3);
        accel_t zi2 = (accel_t)(accz * 2 / 3 + accz_last / 3);
        main_algorithm_trust_ppg4(dm1, pi1, xi1, yi1, zi1);
        main_algorithm_trust_ppg4(dm2, pi2, xi2, yi2, zi2);
        final_HR = main_algorithm_trust_ppg4(delta_ms, ppg, accx, accy, accz);
    }

    ppg_last  = ppg;
    accx_last = accx;
    accy_last = accy;
    accz_last = accz;

    return final_HR;
}

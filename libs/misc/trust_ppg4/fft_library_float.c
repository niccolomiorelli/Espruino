#include "fft_library_float.h"
#include <math.h>
#include <stdbool.h>

#ifndef M_PIF
#define M_PIF 3.14159265358979323846f
#endif

/* Circular buffer index helper */
int buffer_index_plus_fftLib(int buffer_next_i, int plus, int max)
{
    return (buffer_next_i + plus) % max;
}

/* Log base 2 (integer) */
int my_log2_fftLib_float(int N)
{
    int k = N, i = 0;
    while (k) { k >>= 1; i++; }
    return i - 1;
}

/* Bit-reversal of n for log2(N)-bit word */
int reverse_fftLib_float(int N, int n)
{
    int log2N = my_log2_fftLib_float(N);
    int j, p = 0;
    for (j = 1; j <= log2N; j++) {
        if (n & (1 << (log2N - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

/* In-place bit-reversal permutation */
static void bitrev_inplace(complex_number_float *f, int N)
{
    for (int i = 0; i < N; ++i) {
        int j = reverse_fftLib_float(N, i);
        if (j > i) {
            complex_number_float t = f[i];
            f[i] = f[j];
            f[j] = t;
        }
    }
}

/* Pre-computed twiddle factors for N_PAD-point FFT */
static complex_number_float W_half[TP4_N_PAD / 2];
static bool tw_ready = false;

static void ensure_twiddles(void)
{
    if (tw_ready) return;
    for (int k = 0; k < TP4_N_PAD / 2; ++k) {
        float ang = -2.0f * M_PIF * (float)k / (float)TP4_N_PAD;
        W_half[k].real = cosf(ang);
        W_half[k].imag = sinf(ang);
    }
    tw_ready = true;
}

/* Iterative Cooley-Tukey DIT FFT (in-place, N_PAD points) */
static void transform_fftLib_float(complex_number_float *f)
{
    ensure_twiddles();
    bitrev_inplace(f, TP4_N_PAD);

    int step = 1;
    while (step < TP4_N_PAD) {
        int jump   = step << 1;
        int stride = TP4_N_PAD / jump;
        for (int i = 0; i < TP4_N_PAD; i += jump) {
            for (int j = 0; j < step; j++) {
                int i1 = i + j;
                int i2 = i + j + step;
                int W_index = j * stride;

                float c  = W_half[W_index].real;
                float s  = W_half[W_index].imag;

                float tr = c * f[i2].real - s * f[i2].imag;
                float ti = c * f[i2].imag + s * f[i2].real;

                float ur = f[i1].real;
                float ui = f[i1].imag;

                f[i2].real = ur - tr;
                f[i2].imag = ui - ti;
                f[i1].real = ur + tr;
                f[i1].imag = ui + ti;
            }
        }
        step <<= 1;
    }
}

/* Public API */
void FFT_fftLib_float(complex_number_float *f, float d)
{
    transform_fftLib_float(f);
    if (d != 1.0f) {
        for (int i = 0; i < TP4_N_PAD; i++) {
            f[i].real *= d;
            f[i].imag *= d;
        }
    }
}

void apply_hann_window_fftLib_float(const float *input, float *windowed_output, int len)
{
    for (int n = 0; n < len; n++) {
        float hann = 0.5f * (1.0f - cosf(2.0f * M_PIF * (float)n / (float)(len - 1)));
        windowed_output[n] = input[n] * hann;
    }
}

void zero_pad_fftLib_float(const complex_number_float *in, complex_number_float *out, int N)
{
    int i = 0;
    for (; i < N && i < TP4_N_PAD; ++i) out[i] = in[i];
    for (; i < TP4_N_PAD; ++i) { out[i].real = 0.0f; out[i].imag = 0.0f; }
}

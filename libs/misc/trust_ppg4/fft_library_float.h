#ifndef FFT_LIBRARY_FLOAT_H
#define FFT_LIBRARY_FLOAT_H

#include <stdint.h>

#ifndef TP4_WINDOW_LEN
#define TP4_WINDOW_LEN 128
#endif

#ifndef TP4_N_PAD
#define TP4_N_PAD 1024
#endif

typedef struct {
    float real;
    float imag;
} complex_number_float;

/* Circular buffer helper: returns (buffer_next_i + plus) % max */
int buffer_index_plus_fftLib(int buffer_next_i, int plus, int max);

int my_log2_fftLib_float(int N);
int reverse_fftLib_float(int N, int n);

/* In-place FFT of N_PAD points, scaled by d */
void FFT_fftLib_float(complex_number_float *f, float d);

/* Apply Hann window to 'input' of length 'len', write result to 'windowed_output' */
void apply_hann_window_fftLib_float(const float *input, float *windowed_output, int len);

/* Zero-pad 'in' (N points) into 'out' (N_PAD points) */
void zero_pad_fftLib_float(const complex_number_float *in, complex_number_float *out, int N);

#endif

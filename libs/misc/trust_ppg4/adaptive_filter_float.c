#include "adaptive_filter_float.h"

void AdaptFilter_multi_init_float(AdaptFilter_multi_float *f)
{
    int i, j;
    for (i = 0; i < NLMS_ORDER; ++i) {
        for (j = 0; j < NLMS_N_INPUT; j++) {
            f->x[j][i] = 0.0f;
            f->weights[j][i] = 0.0f;
        }
    }
    f->d = 0.0f;
    f->last_index = 0;
}

void AdaptFilter_multi_put_float(AdaptFilter_multi_float *f,
                                  float x_i0, float x_i1, float x_i2,
                                  float d_i)
{
    f->d = d_i;
    /* Shift buffer */
    for (int i = NLMS_ORDER - 1; i > 0; i--) {
        for (int j = 0; j < NLMS_N_INPUT; j++) {
            f->x[j][i] = f->x[j][i - 1];
        }
    }
    f->x[0][0] = x_i0;
    f->x[1][0] = x_i1;
    f->x[2][0] = x_i2;
}

float AdaptFilter_multi_get_float(AdaptFilter_multi_float *f, float mu)
{
    float acc = 0.0f;
    for (int i = 0; i < NLMS_ORDER; i++) {
        for (int j = 0; j < NLMS_N_INPUT; j++) {
            acc += f->weights[j][i] * f->x[j][i];
        }
    }

    float e = f->d - acc;

    float epsilon = 1e-8f;
    float norm = epsilon;
    for (int i = 0; i < NLMS_ORDER; ++i) {
        for (int j = 0; j < NLMS_N_INPUT; j++) {
            norm += f->x[j][i] * f->x[j][i];
        }
    }

    /* NLMS weight update */
    for (int i = 0; i < NLMS_ORDER; ++i) {
        for (int j = 0; j < NLMS_N_INPUT; j++) {
            f->weights[j][i] += (mu / norm) * e * f->x[j][i];
        }
    }

    return e;
}

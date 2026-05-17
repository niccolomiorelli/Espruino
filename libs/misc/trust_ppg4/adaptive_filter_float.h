#ifndef ADAPTIVEFILTER_FLOAT_H_
#define ADAPTIVEFILTER_FLOAT_H_

#include "trust_ppg4_types.h"

#define NLMS_ORDER   75
#define NLMS_N_INPUT  3

/* Single-stage multi-input NLMS adaptive filter */
typedef struct AdaptFilter_multi_float
{
  float x[NLMS_N_INPUT][NLMS_ORDER];
  float weights[NLMS_N_INPUT][NLMS_ORDER];
  float d;
  unsigned int last_index;
} AdaptFilter_multi_float;

void AdaptFilter_multi_init_float(AdaptFilter_multi_float *f);
void AdaptFilter_multi_put_float(AdaptFilter_multi_float *f,
                                  float x_i0, float x_i1, float x_i2,
                                  float d_i);
float AdaptFilter_multi_get_float(AdaptFilter_multi_float *f, float mu);

#endif

#ifndef BPFILTER_H_
#define BPFILTER_H_

#include "trust_ppg4_types.h"

#define BPFILTER_TAP_NUM 35

typedef struct BPFilter
{
  ppg_t history[BPFILTER_TAP_NUM];
  unsigned int last_index;
} BPFilter;

void BPFilter_init(BPFilter *f);
void BPFilter_put(BPFilter *f, ppg_t input);
int BPFilter_get(BPFilter *f);

#endif

/* Rolling statistics using the Welford online algorithm */

#ifndef ROLLING_STATS_FLOAT_H
#define ROLLING_STATS_FLOAT_H

#include <math.h>

#define ROLLING_STATS_WINDOW_SIZE 128

typedef struct Stats_float
{
    unsigned long n;
    int index;
    float old_mean;
    float mean;
    float run_var;
    float run_var_old;
    float window[ROLLING_STATS_WINDOW_SIZE];
    int full;
} Stats_float;

void rolling_stats_reset_float(Stats_float *stats);
void rolling_stats_addValue_float(float x, Stats_float *stats);
float rolling_stats_get_mean_float(Stats_float *stats);
float rolling_stats_get_variance_float(Stats_float *stats);
float rolling_stats_get_standard_deviation_float(Stats_float *stats);

#endif

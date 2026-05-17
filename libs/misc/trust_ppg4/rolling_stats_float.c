#include "rolling_stats_float.h"
#include <math.h>

void rolling_stats_reset_float(Stats_float *stats)
{
    stats->n = 0;
    stats->index = 0;
    stats->old_mean = 0.0f;
    stats->mean = 0.0f;
    stats->run_var = 0.0f;
    stats->run_var_old = 0.0f;
    for (int i = 0; i < ROLLING_STATS_WINDOW_SIZE; i++) {
        stats->window[i] = 0.0f;
    }
    stats->full = 0;
}

void rolling_stats_addValue_float(float x, Stats_float *stats)
{
    int i = stats->index;
    float old_value = stats->window[i];

    stats->window[i] = x;
    stats->index = (i + 1) % ROLLING_STATS_WINDOW_SIZE;

    if (!stats->full) {
        stats->n += 1;
        float delta = x - stats->mean;
        stats->mean += delta / (float)stats->n;
        stats->run_var += delta * (x - stats->mean);
        if (stats->n == ROLLING_STATS_WINDOW_SIZE) {
            stats->full = 1;
        }
    } else {
        stats->old_mean = stats->mean;
        stats->mean += (x - old_value) / (float)ROLLING_STATS_WINDOW_SIZE;
        stats->run_var += (x + old_value - stats->old_mean - stats->mean) * (x - old_value);
        if (stats->run_var < 0.0f) {
            stats->run_var = stats->run_var_old;
        }
        stats->run_var_old = stats->run_var;
    }
}

float rolling_stats_get_mean_float(Stats_float *stats)
{
    return stats->n ? stats->mean : 0.0f;
}

float rolling_stats_get_variance_float(Stats_float *stats)
{
    int denom = stats->full ? ROLLING_STATS_WINDOW_SIZE : (int)stats->n;
    return denom > 1 ? (stats->run_var / (float)(denom - 1)) : 0.0f;
}

float rolling_stats_get_standard_deviation_float(Stats_float *stats)
{
    float variance = rolling_stats_get_variance_float(stats);
    return (float)sqrtf(variance);
}

/* Bandpass Filter
 * Filter designed with scipy.signal firwin()
 * Sampling frequency: 25 Hz, Cut frequencies: (0.5, 4) Hz
 * Order: 35, Window: Hamming, Fixed point precision: 15 bits
 */

#include "bandpass_filter.h"
#include <limits.h>

static int16_t filter_taps[BPFILTER_TAP_NUM] = {
    -89, -72, -28,  2, -82, -336, -626, -661, -308, 104, -85, -1202,
 -2609, -2858, -782, 3298, 7401, 9126, 7401, 3298, -782, -2858, -2609, -1202,
   -85, 104, -308, -661, -626, -336,  -82, 2, -28, -72, -89
};

void BPFilter_init(BPFilter *f)
{
    int i;
    for (i = 0; i < BPFILTER_TAP_NUM; ++i)
        f->history[i] = 0;
    f->last_index = 0;
}

void BPFilter_put(BPFilter *f, ppg_t input)
{
    f->history[f->last_index++] = input;
    if (f->last_index == BPFILTER_TAP_NUM)
        f->last_index = 0;
}

int BPFilter_get(BPFilter *f)
{
    long long acc = 0;
    int index = f->last_index, i;
    for (i = 0; i < BPFILTER_TAP_NUM; ++i)
    {
        index = index != 0 ? index - 1 : BPFILTER_TAP_NUM - 1;
        acc += (long long)f->history[index] * filter_taps[i];
    }
    int temp = (int)(acc >> 15);
    if (temp > INT_MAX) temp = INT_MAX;
    if (temp < INT_MIN) temp = INT_MIN;
    return temp;
}

/**
 * @file bench_util.c
 * @brief Timing and statistics helpers shared by the benchmark harness.
 */

#include "bench.h"

#include <time.h>
#include <stdlib.h>
#include <math.h>

/**
 * Monotonic wall-clock time in seconds.
 *
 * Uses CLOCK_MONOTONIC so the measurement is unaffected by wall-clock
 * adjustments; the generation routines are multi-threaded, so wall time is the
 * quantity of interest (the real latency to grow one CFF into the next).
 *
 * @return Seconds since an unspecified fixed point, as a double.
 */
double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/**
 * Comparator for qsort of doubles (ascending).
 */
static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/**
 * Computes summary statistics (mean, min, max, sample std-dev, median) over a
 * set of timing samples.
 *
 * @param x Array of samples.
 * @param n Number of samples.
 * @return Filled Stats; a zeroed struct when n <= 0.
 */
Stats stats_compute(const double *x, int n)
{
    Stats s = (Stats){0};
    s.n = n;
    if (n <= 0) return s;

    double sum = 0.0, mn = x[0], mx = x[0];
    for (int i = 0; i < n; i++) {
        sum += x[i];
        if (x[i] < mn) mn = x[i];
        if (x[i] > mx) mx = x[i];
    }
    s.mean = sum / n;
    s.min  = mn;
    s.max  = mx;

    if (n >= 2) {
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            double dlt = x[i] - s.mean;
            acc += dlt * dlt;
        }
        s.std = sqrt(acc / (n - 1));   /* sample standard deviation */
    }

    double *tmp = malloc((size_t)n * sizeof(double));
    if (tmp) {
        for (int i = 0; i < n; i++) tmp[i] = x[i];
        qsort(tmp, (size_t)n, sizeof(double), cmp_double);
        s.median = (n & 1) ? tmp[n / 2]
                           : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
        free(tmp);
    }
    return s;
}

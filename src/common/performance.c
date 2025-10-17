#include "performance.h"

#include <stddef.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void ns_perf_now(struct timespec *ts) {
#if defined(CLOCK_MONOTONIC_RAW)
    clock_gettime(CLOCK_MONOTONIC_RAW, ts);
#elif defined(CLOCK_MONOTONIC)
    clock_gettime(CLOCK_MONOTONIC, ts);
#else
    clock_gettime(CLOCK_REALTIME, ts);
#endif
}

void ns_perf_timer_start(NSPERFTimer *timer) {
    if (!timer) return;
    ns_perf_now(&timer->start);
    timer->running = 1;
}

double ns_perf_timer_stop(NSPERFTimer *timer) {
    if (!timer) return 0.0;
    ns_perf_now(&timer->stop);
    timer->running = 0;
    return ns_perf_elapsed_seconds(&timer->start, &timer->stop);
}

double ns_perf_elapsed_seconds(const struct timespec *start, const struct timespec *stop) {
    if (!start || !stop) return 0.0;
    const double sec = (double)(stop->tv_sec - start->tv_sec);
    const double nsec = (double)(stop->tv_nsec - start->tv_nsec) / 1.0e9;
    return sec + nsec;
}

int ns_perf_active_threads(void) {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

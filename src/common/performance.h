#ifndef NS_COMMON_PERFORMANCE_H
#define NS_COMMON_PERFORMANCE_H

#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    struct timespec start;
    struct timespec stop;
    int running;
} NSPERFTimer;

void ns_perf_timer_start(NSPERFTimer *timer);
double ns_perf_timer_stop(NSPERFTimer *timer);
double ns_perf_elapsed_seconds(const struct timespec *start, const struct timespec *stop);
int ns_perf_active_threads(void);

#ifdef __cplusplus
}
#endif

#endif /* NS_COMMON_PERFORMANCE_H */

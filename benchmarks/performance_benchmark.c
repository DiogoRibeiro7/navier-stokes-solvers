#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/common/performance.h"
#include "../src/finite_difference/ns_fd_solver.h"
#include "../src/spectral/ns_spectral_solver.h"

typedef struct {
    int nx;
    int ny;
    int steps;
    double Re;
} bench_case;

static void run_fd_case(const bench_case *cfg, FILE *out) {
    NSFiniteDiffData *data = ns_fd_allocate(cfg->nx, cfg->ny, 1.0, 1.0, cfg->Re);
    if (!data) {
        fprintf(stderr, "[bench] failed to allocate FD solver (%dx%d)\n", cfg->nx, cfg->ny);
        return;
    }
    ns_fd_initialize_lid_cavity(data);

    ns_perf_counter timer;
    const size_t bytes_per_step =
        (size_t)cfg->nx * (size_t)cfg->ny * sizeof(double) * 6;
    ns_perf_counter_start(&timer, bytes_per_step * (size_t)cfg->steps, 0);

    for (int step = 0; step < cfg->steps; ++step) {
        if (ns_fd_advance_timestep(data) == SOLVER_DIVERGED) {
            fprintf(stderr, "[bench] FD solver diverged at step %d\n", step);
            break;
        }
    }

    ns_perf_counter_stop(&timer);
    const double wall = ns_perf_seconds(&timer.start, &timer.end);
    fprintf(out,
            "fd,%d,%d,%d,%d,%.6f,%.3f,%.2f\n",
            cfg->nx,
            cfg->ny,
            cfg->steps,
            ns_perf_active_threads(),
            wall,
            ns_perf_counter_bandwidth_gbps(&timer),
            ns_perf_resident_memory() / (1024.0 * 1024.0));

    ns_fd_free(data);
}

static void run_spectral_case(const bench_case *cfg, FILE *out) {
    NSSpectralData *data =
        ns_spectral_allocate(cfg->nx, cfg->ny, 2.0 * M_PI, 2.0 * M_PI, cfg->Re);
    if (!data) {
        fprintf(stderr, "[bench] failed to allocate spectral solver (%dx%d)\n",
                cfg->nx, cfg->ny);
        return;
    }

    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    ns_spectral_compute_derivatives(data);

    ns_perf_counter timer;
    const size_t bytes_per_step =
        (size_t)cfg->nx * (size_t)cfg->ny * sizeof(double) * 8;
    ns_perf_counter_start(&timer, bytes_per_step * (size_t)cfg->steps, 0);

    for (int step = 0; step < cfg->steps; ++step) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
    }

    ns_perf_counter_stop(&timer);
    const double wall = ns_perf_seconds(&timer.start, &timer.end);
    fprintf(out,
            "spectral,%d,%d,%d,%d,%.6f,%.3f,%.2f\n",
            cfg->nx,
            cfg->ny,
            cfg->steps,
            ns_perf_active_threads(),
            wall,
            ns_perf_counter_bandwidth_gbps(&timer),
            ns_perf_resident_memory() / (1024.0 * 1024.0));

    ns_spectral_free(data);
}

int main(void) {
    ns_perf_configure_threading("Performance benchmark");
    printf("solver,nx,ny,steps,threads,wall_time_s,bandwidth_gbps,resident_mem_mib\n");

    const bench_case fd_cases[] = {
        {64, 64, 200, 100.0},
        {96, 96, 150, 100.0},
        {128, 128, 100, 400.0},
    };

    const bench_case spec_cases[] = {
        {64, 64, 400, 800.0},
        {128, 128, 250, 1600.0},
        {192, 192, 150, 1600.0},
    };

    for (size_t i = 0; i < sizeof(fd_cases) / sizeof(fd_cases[0]); ++i) {
        run_fd_case(&fd_cases[i], stdout);
    }
    for (size_t i = 0; i < sizeof(spec_cases) / sizeof(spec_cases[0]); ++i) {
        run_spectral_case(&spec_cases[i], stdout);
    }

    return 0;
}

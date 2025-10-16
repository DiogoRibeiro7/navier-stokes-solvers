#include "ns_fd_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Compute solver statistics
void ns_fd_compute_stats(NSFiniteDiffData *data, SolverStats *stats) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    
    double energy_sum = 0.0;
    double max_vel = 0.0;
    double max_div = 0.0;

    const double *NS_RESTRICT u = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->u, 64);
    const double *NS_RESTRICT v = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->v, 64);

#pragma omp parallel for collapse(2) reduction(+:energy_sum) reduction(max:max_vel, max_div)
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            const int ij = IDX2D(i, j, nx);
            const double vel_mag = sqrt(u[ij] * u[ij] + v[ij] * v[ij]);
            energy_sum += 0.5 * vel_mag * vel_mag;
            if (vel_mag > max_vel) {
                max_vel = vel_mag;
            }

            const double du_dx = (u[IDX2D(i, j + 1, nx)] - u[IDX2D(i, j - 1, nx)]) / (2 * dx);
            const double dv_dy = (v[IDX2D(i + 1, j, nx)] - v[IDX2D(i - 1, j, nx)]) / (2 * dy);
            const double div = fabs(du_dx + dv_dy);
            if (div > max_div) {
                max_div = div;
            }
        }
    }
    
    stats->energy = energy_sum / (nx * ny);
    stats->max_velocity = max_vel;
    stats->max_divergence = max_div;
    stats->residual = ns_fd_compute_residual(data);
}

// Analyze convergence
void ns_fd_analyze_convergence(NSFiniteDiffData *data, int step) {
    if (step % 10 != 0) return;
    
    SolverStats stats;
    ns_fd_compute_stats(data, &stats);
    
    printf("Step %d: t=%.4f, E=%.6e, max_vel=%.4f, max_div=%.2e, res=%.2e, dt=%.2e\n",
           step, data->t, stats.energy, stats.max_velocity, 
           stats.max_divergence, stats.residual, data->dt);
}

// Output solution to file
void ns_fd_output_solution(NSFiniteDiffData *data, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return;
    
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    
    fprintf(fp, "# x y u v p\n");
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            int ij = IDX2D(i, j, nx);
            double x = j * dx;
            double y = i * dy;
            fprintf(fp, "%.6f %.6f %.6e %.6e %.6e\n",
                   x, y, data->u[ij], data->v[ij], data->p[ij]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

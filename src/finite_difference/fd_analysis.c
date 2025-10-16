#include "ns_fd_solver.h"
#include <time.h>

// Compute solver statistics
void ns_fd_compute_stats(NSFiniteDiffData *data, SolverStats *stats) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    
    stats->energy = 0.0;
    stats->max_velocity = 0.0;
    stats->max_divergence = 0.0;
    
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            double vel_mag = sqrt(data->u[ij]*data->u[ij] + data->v[ij]*data->v[ij]);
            stats->energy += 0.5 * vel_mag * vel_mag;
            stats->max_velocity = fmax(stats->max_velocity, vel_mag);
            
            double du_dx = (data->u[IDX2D(i, j+1, nx)] - data->u[IDX2D(i, j-1, nx)]) / (2*dx);
            double dv_dy = (data->v[IDX2D(i+1, j, nx)] - data->v[IDX2D(i-1, j, nx)]) / (2*dy);
            double div = fabs(du_dx + dv_dy);
            stats->max_divergence = fmax(stats->max_divergence, div);
        }
    }
    
    stats->energy /= (nx * ny);
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

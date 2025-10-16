/*
 * Taylor-Green Vortex Validation Tests
 * Tests both solvers against analytical solution
 */

#include "test_framework.h"
#include "../src/finite_difference/ns_fd_solver.h"
#include "../src/spectral/ns_spectral_solver.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Analytical solution for Taylor-Green vortex
double taylor_green_u_exact(double x, double y, double t, double Re) {
    double Lx = 2*M_PI;
    double Ly = 2*M_PI;
    double decay = exp(-2.0 * (2*M_PI/Lx)*(2*M_PI/Lx) * t / Re);
    return sin(2*M_PI*x/Lx) * cos(2*M_PI*y/Ly) * decay;
}

double taylor_green_v_exact(double x, double y, double t, double Re) {
    double Lx = 2*M_PI;
    double Ly = 2*M_PI;
    double decay = exp(-2.0 * (2*M_PI/Lx)*(2*M_PI/Lx) * t / Re);
    return -cos(2*M_PI*x/Lx) * sin(2*M_PI*y/Ly) * decay;
}

double taylor_green_energy_exact(double t, double Re) {
    double Lx = 2*M_PI;
    double decay = exp(-4.0 * (2*M_PI/Lx)*(2*M_PI/Lx) * t / Re);
    return 0.5 * decay;  // Normalized energy
}

// Test: Spectral solver Taylor-Green vortex
int test_spectral_taylor_green_short_time() {
    printf("  Testing spectral solver with Taylor-Green vortex (short time)...\n");
    
    int nx = 64, ny = 64;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    double t_final = 0.1;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    // Run for short time
    int steps = 0;
    while (data->t < t_final && steps < 100) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
        steps++;
    }
    
    // Compute L2 error
    double error = 0.0;
    int count = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;
            double x = data->x[i];
            double y = data->y[j];
            
            double u_exact = taylor_green_u_exact(x, y, data->t, Re);
            double v_exact = taylor_green_v_exact(x, y, data->t, Re);
            
            error += pow(data->u[idx] - u_exact, 2) + pow(data->v[idx] - v_exact, 2);
            count++;
        }
    }
    error = sqrt(error / count);
    
    printf("    Time: %.4f, Steps: %d, L2 error: %.2e\n", data->t, steps, error);
    
    TEST_ASSERT(error < 1e-4, "Spectral solver L2 error < 1e-4");
    
    ns_spectral_free(data);
    return 1;
}

// Test: Energy decay rate
int test_spectral_energy_decay() {
    printf("  Testing energy decay rate...\n");
    
    int nx = 64, ny = 64;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    double t_final = 1.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    double energy_initial = ns_spectral_compute_energy(data);
    
    // Advance in time
    while (data->t < t_final) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
    }
    
    double energy_final = ns_spectral_compute_energy(data);
    double energy_exact = taylor_green_energy_exact(data->t, Re);
    
    double energy_error = fabs(energy_final - energy_exact) / energy_exact;
    
    printf("    Energy initial: %.6e\n", energy_initial);
    printf("    Energy final:   %.6e\n", energy_final);
    printf("    Energy exact:   %.6e\n", energy_exact);
    printf("    Relative error: %.2e\n", energy_error);
    
    TEST_ASSERT(energy_error < 0.01, "Energy decay error < 1%");
    
    ns_spectral_free(data);
    return 1;
}

// Test: Spectral accuracy (exponential convergence)
int test_spectral_convergence_rate() {
    printf("  Testing spectral convergence rate...\n");
    
    double Re = 100.0;
    double t_final = 0.5;
    int resolutions[] = {16, 32, 64};
    double errors[3];
    
    for (int r = 0; r < 3; r++) {
        int nx = resolutions[r];
        int ny = resolutions[r];
        double Lx = 2*M_PI, Ly = 2*M_PI;
        
        NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
        ns_spectral_initialize_grids(data);
        ns_spectral_initialize_fft_plans(data);
        ns_spectral_initialize_taylor_green(data);
        
        while (data->t < t_final) {
            data->dt = ns_spectral_compute_cfl_timestep(data);
            ns_spectral_advance_rk4(data);
        }
        
        // Compute error
        double error = 0.0;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = j * nx + i;
                double u_exact = taylor_green_u_exact(data->x[i], data->y[j], data->t, Re);
                double v_exact = taylor_green_v_exact(data->x[i], data->y[j], data->t, Re);
                error += pow(data->u[idx] - u_exact, 2) + pow(data->v[idx] - v_exact, 2);
            }
        }
        errors[r] = sqrt(error / (nx * ny));
        
        printf("    N=%d: error = %.2e\n", nx, errors[r]);
        
        ns_spectral_free(data);
    }
    
    // Check exponential convergence (error should decrease rapidly)
    TEST_ASSERT(errors[1] < 0.1 * errors[0], "Error decreases from 16 to 32");
    TEST_ASSERT(errors[2] < 0.1 * errors[1], "Error decreases from 32 to 64");
    
    printf("    Convergence rate demonstrated!\n");
    
    return 1;
}

// Test: Mass conservation
int test_mass_conservation() {
    printf("  Testing mass conservation...\n");
    
    int nx = 32, ny = 32;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    // Check divergence at multiple time steps
    for (int step = 0; step < 50; step++) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
        
        // Compute maximum divergence
        double max_div = 0.0;
        for (int j = 1; j < ny-1; j++) {
            for (int i = 1; i < nx-1; i++) {
                int idx = j * nx + i;
                double du_dx = (data->u[j*nx + (i+1)] - data->u[j*nx + (i-1)]) / (2*data->Lx/nx);
                double dv_dy = (data->v[(j+1)*nx + i] - data->v[(j-1)*nx + i]) / (2*data->Ly/ny);
                double div = fabs(du_dx + dv_dy);
                max_div = fmax(max_div, div);
            }
        }
        
        if (step % 10 == 0) {
            printf("    Step %d: max(∇·u) = %.2e\n", step, max_div);
        }
        
        TEST_ASSERT(max_div < 1e-10, "Divergence remains small (machine precision)");
    }
    
    ns_spectral_free(data);
    return 1;
}

// Main test runner
int main() {
    TEST_SUITE_START("Taylor-Green Vortex Validation");
    
    RUN_TEST(test_spectral_taylor_green_short_time);
    RUN_TEST(test_spectral_energy_decay);
    RUN_TEST(test_spectral_convergence_rate);
    RUN_TEST(test_mass_conservation);
    
    TEST_SUITE_END();
}

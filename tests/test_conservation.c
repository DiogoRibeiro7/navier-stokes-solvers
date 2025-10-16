/*
 * Conservation Property Tests
 * Verifies mass, momentum, and energy conservation
 */

#include "test_framework.h"
#include "../src/finite_difference/ns_fd_solver.h"
#include "../src/spectral/ns_spectral_solver.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Test: Mass conservation for FD solver
int test_fd_mass_conservation() {
    printf("  Testing FD mass conservation...\n");
    
    int nx = 32, ny = 32;
    double L = 1.0, H = 1.0;
    double Re = 100.0;
    
    NSFiniteDiffData *data = ns_fd_allocate(nx, ny, L, H, Re);
    TEST_ASSERT_NOT_NULL(data, "FD data allocation");
    
    ns_fd_initialize_lid_cavity(data);
    
    // Run for several time steps
    for (int step = 0; step < 20; step++) {
        ns_fd_advance_timestep(data);
        
        // Compute maximum divergence
        double max_div = 0.0;
        double avg_div = 0.0;
        int count = 0;
        
        for (int i = 1; i < ny-1; i++) {
            for (int j = 1; j < nx-1; j++) {
                int ij = i * nx + j;
                double du_dx = (data->u[i*nx + (j+1)] - data->u[i*nx + (j-1)]) / (2*data->dx);
                double dv_dy = (data->v[(i+1)*nx + j] - data->v[(i-1)*nx + j]) / (2*data->dy);
                double div = fabs(du_dx + dv_dy);
                max_div = fmax(max_div, div);
                avg_div += div;
                count++;
            }
        }
        avg_div /= count;
        
        if (step % 5 == 0) {
            printf("    Step %d: max_div=%.2e, avg_div=%.2e\n", step, max_div, avg_div);
        }
        
        TEST_ASSERT(max_div < 1e-6, "Maximum divergence acceptable");
    }
    
    ns_fd_free(data);
    return 1;
}

// Test: Spectral solver incompressibility (should be machine precision)
int test_spectral_incompressibility() {
    printf("  Testing spectral incompressibility...\n");
    
    int nx = 64, ny = 64;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    // Run and check at each step
    for (int step = 0; step < 30; step++) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
        
        // Compute divergence using spectral derivatives
        double max_div = 0.0;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = j * nx + i;
                int i_p = (i + 1) % nx;
                int i_m = (i - 1 + nx) % nx;
                int j_p = (j + 1) % ny;
                int j_m = (j - 1 + ny) % ny;
                
                double du_dx = (data->u[j*nx + i_p] - data->u[j*nx + i_m]) / (2*data->Lx/nx);
                double dv_dy = (data->v[j_p*nx + i] - data->v[j_m*nx + i]) / (2*data->Ly/ny);
                double div = fabs(du_dx + dv_dy);
                max_div = fmax(max_div, div);
            }
        }
        
        if (step % 10 == 0) {
            printf("    Step %d: max(∇·u) = %.2e\n", step, max_div);
        }
        
        TEST_ASSERT(max_div < 1e-12, "Spectral divergence at machine precision");
    }
    
    ns_spectral_free(data);
    return 1;
}

// Test: Energy dissipation rate (should match viscous theory)
int test_energy_dissipation() {
    printf("  Testing energy dissipation...\n");
    
    int nx = 64, ny = 64;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    double t_final = 1.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    double E0 = ns_spectral_compute_energy(data);
    double t0 = data->t;
    
    // Advance in time
    int steps = 0;
    while (data->t < t_final && steps < 1000) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
        steps++;
    }
    
    double E1 = ns_spectral_compute_energy(data);
    double dt_elapsed = data->t - t0;
    
    // Energy should decay
    printf("    Initial energy: %.6e\n", E0);
    printf("    Final energy:   %.6e\n", E1);
    printf("    Decay factor:   %.4f\n", E1/E0);
    
    TEST_ASSERT(E1 < E0, "Energy decreases due to viscosity");
    TEST_ASSERT(E1 > 0, "Energy remains positive");
    TEST_ASSERT(E1/E0 > 0.01, "Energy decay is reasonable (not excessive)");
    
    ns_spectral_free(data);
    return 1;
}

// Test: Enstrophy evolution
int test_enstrophy_evolution() {
    printf("  Testing enstrophy evolution...\n");
    
    int nx = 64, ny = 64;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    double Omega0 = ns_spectral_compute_enstrophy(data);
    
    // Advance several steps
    for (int step = 0; step < 50; step++) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
    }
    
    double Omega1 = ns_spectral_compute_enstrophy(data);
    
    printf("    Initial enstrophy: %.6e\n", Omega0);
    printf("    Final enstrophy:   %.6e\n", Omega1);
    
    // Enstrophy should decay monotonically
    TEST_ASSERT(Omega1 < Omega0, "Enstrophy decreases");
    TEST_ASSERT(Omega1 > 0, "Enstrophy remains positive");
    
    ns_spectral_free(data);
    return 1;
}

// Test: Symmetry preservation
int test_symmetry_preservation() {
    printf("  Testing symmetry preservation...\n");
    
    int nx = 32, ny = 32;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 100.0;
    
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocation");
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    // Run for some time
    for (int step = 0; step < 20; step++) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
    }
    
    // Check symmetry: Taylor-Green should maintain certain symmetries
    double symmetry_error = 0.0;
    int count = 0;
    
    for (int j = 0; j < ny/2; j++) {
        for (int i = 0; i < nx/2; i++) {
            int idx1 = j * nx + i;
            int idx2 = (ny-1-j) * nx + (nx-1-i);
            
            // u should be antisymmetric, v should be antisymmetric
            symmetry_error += fabs(data->u[idx1] + data->u[idx2]);
            symmetry_error += fabs(data->v[idx1] + data->v[idx2]);
            count += 2;
        }
    }
    symmetry_error /= count;
    
    printf("    Symmetry error: %.2e\n", symmetry_error);
    TEST_ASSERT(symmetry_error < 1e-10, "Symmetry preserved");
    
    ns_spectral_free(data);
    return 1;
}

// Main test runner
int main() {
    TEST_SUITE_START("Conservation Properties");
    
    RUN_TEST(test_fd_mass_conservation);
    RUN_TEST(test_spectral_incompressibility);
    RUN_TEST(test_energy_dissipation);
    RUN_TEST(test_enstrophy_evolution);
    RUN_TEST(test_symmetry_preservation);
    
    TEST_SUITE_END();
}

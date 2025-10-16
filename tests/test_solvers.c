/*
 * Basic Solver Functionality Tests
 * Tests initialization, memory, and basic operations
 */

#include "test_framework.h"
#include "../src/finite_difference/ns_fd_solver.h"
#include "../src/spectral/ns_spectral_solver.h"

// FD Solver Tests
int test_fd_allocation() {
    printf("  Testing FD memory allocation...\n");
    
    NSFiniteDiffData *data = ns_fd_allocate(32, 32, 1.0, 1.0, 100.0);
    TEST_ASSERT_NOT_NULL(data, "FD data structure allocated");
    TEST_ASSERT_NOT_NULL(data->u, "Velocity u allocated");
    TEST_ASSERT_NOT_NULL(data->v, "Velocity v allocated");
    TEST_ASSERT_NOT_NULL(data->p, "Pressure allocated");
    TEST_ASSERT_EQUAL(data->nx, 32, "Grid size nx correct");
    TEST_ASSERT_EQUAL(data->ny, 32, "Grid size ny correct");
    
    ns_fd_free(data);
    return 1;
}

int test_fd_initialization() {
    printf("  Testing FD initialization...\n");
    
    NSFiniteDiffData *data = ns_fd_allocate(32, 32, 1.0, 1.0, 100.0);
    TEST_ASSERT_NOT_NULL(data, "FD data allocated");
    
    ns_fd_initialize_lid_cavity(data);
    
    // Check lid velocity
    int top_row = (data->ny - 1) * data->nx;
    double lid_velocity = data->u[top_row + data->nx/2];
    TEST_ASSERT_DOUBLE_EQUAL(lid_velocity, 1.0, 1e-10, "Lid velocity set correctly");
    
    // Check bottom wall (should be zero)
    double bottom_velocity = data->u[data->nx/2];
    TEST_ASSERT_DOUBLE_EQUAL(bottom_velocity, 0.0, 1e-10, "Bottom wall is no-slip");
    
    ns_fd_free(data);
    return 1;
}

int test_fd_boundary_conditions() {
    printf("  Testing FD boundary conditions...\n");
    
    NSFiniteDiffData *data = ns_fd_allocate(32, 32, 1.0, 1.0, 100.0);
    ns_fd_initialize_lid_cavity(data);
    
    // Modify interior and reapply BCs
    for (int i = 0; i < data->ny; i++) {
        for (int j = 0; j < data->nx; j++) {
            int ij = i * data->nx + j;
            data->u[ij] = 0.5;
            data->v[ij] = 0.3;
        }
    }
    
    ns_fd_apply_boundary_conditions(data);
    
    // Check all walls
    for (int j = 0; j < data->nx; j++) {
        TEST_ASSERT_DOUBLE_EQUAL(data->v[j], 0.0, 1e-10, "Bottom v-velocity zero");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[(data->ny-1)*data->nx + j], 0.0, 1e-10, "Top v-velocity zero");
    }
    
    ns_fd_free(data);
    return 1;
}

int test_fd_timestep() {
    printf("  Testing FD adaptive timestep...\n");
    
    NSFiniteDiffData *data = ns_fd_allocate(32, 32, 1.0, 1.0, 100.0);
    ns_fd_initialize_lid_cavity(data);
    
    double dt = ns_fd_adaptive_timestep(data);
    
    TEST_ASSERT(dt > 0, "Timestep is positive");
    TEST_ASSERT(dt < 1.0, "Timestep is reasonable");
    
    printf("    Computed dt: %.4e\n", dt);
    
    ns_fd_free(data);
    return 1;
}

int test_fd_residual() {
    printf("  Testing FD residual computation...\n");
    
    NSFiniteDiffData *data = ns_fd_allocate(32, 32, 1.0, 1.0, 100.0);
    ns_fd_initialize_lid_cavity(data);
    
    // Store initial values
    memcpy(data->u_old, data->u, data->nx * data->ny * sizeof(double));
    memcpy(data->v_old, data->v, data->nx * data->ny * sizeof(double));
    memcpy(data->p_old, data->p, data->nx * data->ny * sizeof(double));
    
    double residual = ns_fd_compute_residual(data);
    
    TEST_ASSERT(residual >= 0, "Residual is non-negative");
    printf("    Initial residual: %.2e\n", residual);
    
    ns_fd_free(data);
    return 1;
}

// Spectral Solver Tests
int test_spectral_allocation() {
    printf("  Testing spectral memory allocation...\n");
    
    NSSpectralData *data = ns_spectral_allocate(64, 64, 2*M_PI, 2*M_PI, 100.0);
    TEST_ASSERT_NOT_NULL(data, "Spectral data structure allocated");
    TEST_ASSERT_NOT_NULL(data->u, "Velocity u allocated");
    TEST_ASSERT_NOT_NULL(data->v, "Velocity v allocated");
    TEST_ASSERT_NOT_NULL(data->omega, "Vorticity allocated");
    TEST_ASSERT_NOT_NULL(data->u_hat, "Spectral u_hat allocated");
    TEST_ASSERT_NOT_NULL(data->omega_hat, "Spectral omega_hat allocated");
    TEST_ASSERT_EQUAL(data->nx, 64, "Grid size nx correct");
    TEST_ASSERT_EQUAL(data->nkx, 33, "Spectral size nkx correct (nx/2+1)");
    
    ns_spectral_free(data);
    return 1;
}

int test_spectral_grid_initialization() {
    printf("  Testing spectral grid initialization...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocated");
    
    ns_spectral_initialize_grids(data);
    
    // Check wavenumber arrays
    TEST_ASSERT_DOUBLE_EQUAL(data->kx[0], 0.0, 1e-10, "Zero wavenumber at k=0");
    TEST_ASSERT_DOUBLE_EQUAL(data->kx[1], 1.0, 1e-10, "First wavenumber is 1");
    
    // Check k^2 array
    double k2_00 = data->k2[0];
    TEST_ASSERT_DOUBLE_EQUAL(k2_00, 0.0, 1e-10, "k^2 is zero at origin");
    
    // Check dealiasing mask
    TEST_ASSERT_EQUAL(data->dealias_mask[0], 1, "Origin not dealiased");
    
    ns_spectral_free(data);
    return 1;
}

int test_spectral_fft_transforms() {
    printf("  Testing FFT forward/backward transforms...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    
    // Create a known signal
    for (int j = 0; j < data->ny; j++) {
        for (int i = 0; i < data->nx; i++) {
            int idx = j * data->nx + i;
            data->omega[idx] = sin(data->x[i]) * cos(data->y[j]);
        }
    }
    
    // Store original
    double *original = malloc(data->nx * data->ny * sizeof(double));
    memcpy(original, data->omega, data->nx * data->ny * sizeof(double));
    
    // Forward transform
    ns_spectral_physical_to_spectral(data, data->omega, data->omega_hat);
    
    // Backward transform
    ns_spectral_spectral_to_physical(data, data->omega_hat, data->omega);
    
    // Check roundtrip accuracy
    double max_error = 0.0;
    for (int k = 0; k < data->nx * data->ny; k++) {
        double error = fabs(data->omega[k] - original[k]);
        max_error = fmax(max_error, error);
    }
    
    printf("    FFT roundtrip error: %.2e\n", max_error);
    TEST_ASSERT(max_error < 1e-10, "FFT roundtrip accurate");
    
    free(original);
    ns_spectral_free(data);
    return 1;
}

int test_spectral_dealiasing() {
    printf("  Testing dealiasing...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    
    // Create field in spectral space
    for (int k = 0; k < data->nkx * data->nky; k++) {
        data->omega_hat[k] = 1.0 + 0.0*I;
    }
    
    // Apply dealiasing
    ns_spectral_apply_dealiasing(data, data->omega_hat);
    
    // Check that high wavenumbers are zeroed
    int high_k_count = 0;
    int zeroed_count = 0;
    for (int k = 0; k < data->nkx * data->nky; k++) {
        if (!data->dealias_mask[k]) {
            high_k_count++;
            if (cabs(data->omega_hat[k]) < 1e-15) {
                zeroed_count++;
            }
        }
    }
    
    printf("    High-k modes: %d, Zeroed: %d\n", high_k_count, zeroed_count);
    TEST_ASSERT_EQUAL(high_k_count, zeroed_count, "All high-k modes zeroed");
    
    ns_spectral_free(data);
    return 1;
}

int test_spectral_velocity_from_vorticity() {
    printf("  Testing velocity from vorticity...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    
    // Set uniform vorticity in spectral space (k=0 mode only)
    for (int k = 0; k < data->nkx * data->nky; k++) {
        data->omega_hat[k] = 0.0 + 0.0*I;
    }
    data->omega_hat[0] = 1.0 + 0.0*I;  // Uniform vorticity
    
    ns_spectral_compute_velocity_from_vorticity(data);
    
    // For uniform vorticity, velocity should be zero (no circulation)
    ns_spectral_spectral_to_physical(data, data->u_hat, data->u);
    ns_spectral_spectral_to_physical(data, data->v_hat, data->v);
    
    double max_vel = 0.0;
    for (int k = 0; k < data->nx * data->ny; k++) {
        max_vel = fmax(max_vel, fabs(data->u[k]));
        max_vel = fmax(max_vel, fabs(data->v[k]));
    }
    
    printf("    Max velocity from uniform vorticity: %.2e\n", max_vel);
    TEST_ASSERT(max_vel < 1e-10, "Uniform vorticity gives zero velocity");
    
    ns_spectral_free(data);
    return 1;
}

int test_spectral_energy_computation() {
    printf("  Testing energy computation...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    double energy = ns_spectral_compute_energy(data);
    
    TEST_ASSERT(energy > 0, "Energy is positive");
    TEST_ASSERT(energy < 1.0, "Energy is bounded");
    
    printf("    Initial energy: %.6e\n", energy);
    
    ns_spectral_free(data);
    return 1;
}

int test_spectral_enstrophy_computation() {
    printf("  Testing enstrophy computation...\n");
    
    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    double enstrophy = ns_spectral_compute_enstrophy(data);
    
    TEST_ASSERT(enstrophy > 0, "Enstrophy is positive");
    
    printf("    Initial enstrophy: %.6e\n", enstrophy);
    
    ns_spectral_free(data);
    return 1;
}

// Main test runner
int main() {
    TEST_SUITE_START("Solver Functionality Tests");
    
    // FD tests
    RUN_TEST(test_fd_allocation);
    RUN_TEST(test_fd_initialization);
    RUN_TEST(test_fd_boundary_conditions);
    RUN_TEST(test_fd_timestep);
    RUN_TEST(test_fd_residual);
    
    // Spectral tests
    RUN_TEST(test_spectral_allocation);
    RUN_TEST(test_spectral_grid_initialization);
    RUN_TEST(test_spectral_fft_transforms);
    RUN_TEST(test_spectral_dealiasing);
    RUN_TEST(test_spectral_velocity_from_vorticity);
    RUN_TEST(test_spectral_energy_computation);
    RUN_TEST(test_spectral_enstrophy_computation);
    
    TEST_SUITE_END();
}

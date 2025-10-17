/*
 * Basic Solver Functionality Tests
 * Tests initialization, memory, and basic operations
 */

#include "test_framework.h"
#include "../include/ns_config.h"
#include "../src/common/forcing.h"
#include "../src/finite_difference/ns_fd_solver.h"
#include "../src/finite_difference/multigrid.h"
#include "../src/spectral/ns_spectral_solver.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double u;
    double v;
    double p;
} CustomBCData;

static void custom_bc_callback(double x, double y, double t,
                               double *u, double *v, double *p,
                               void *user_data) {
    (void)x;
    (void)y;
    (void)t;
    const CustomBCData *params = (const CustomBCData *)user_data;
    if (!params) {
        return;
    }
    *u = params->u;
    *v = params->v;
    *p = params->p;
}

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
    TEST_ASSERT_NOT_NULL(data, "FD data allocated for BC tests");

    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_NO_SLIP, 1.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_RIGHT, BC_NO_SLIP, 0.0, 0.0, NAN);

    ns_fd_initialize_lid_cavity(data);

    for (int i = 0; i < data->ny; ++i) {
        for (int j = 0; j < data->nx; ++j) {
            const int ij = i * data->nx + j;
            data->u[ij] = 0.5;
            data->v[ij] = 0.3;
        }
    }

    ns_fd_apply_boundary_conditions(data);

    for (int j = 0; j < data->nx; ++j) {
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(0, j, data->nx)], 0.0, 1e-10, "Bottom v-velocity zero for no-slip");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(data->ny - 1, j, data->nx)], 0.0, 1e-10, "Top v-velocity zero for no-slip");
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(data->ny - 1, j, data->nx)], 1.0, 1e-10, "Top u-velocity matches lid speed");
    }

    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_FREE_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_FREE_SLIP, 0.0, 0.0, NAN);

    memset(data->u, 0, (size_t)data->nx * data->ny * sizeof(double));
    memset(data->v, 0, (size_t)data->nx * data->ny * sizeof(double));

    for (int j = 0; j < data->nx; ++j) {
        data->u[IDX2D(1, j, data->nx)] = 0.2 + 0.01 * j;
        data->v[IDX2D(1, j, data->nx)] = 0.05;
        data->u[IDX2D(data->ny - 2, j, data->nx)] = -0.3 + 0.02 * j;
        data->v[IDX2D(data->ny - 2, j, data->nx)] = -0.04;
    }

    ns_fd_apply_boundary_conditions(data);

    for (int j = 0; j < data->nx; ++j) {
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(0, j, data->nx)], data->u[IDX2D(1, j, data->nx)], 1e-10, "Free-slip bottom copies interior u");
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(data->ny - 1, j, data->nx)], data->u[IDX2D(data->ny - 2, j, data->nx)], 1e-10, "Free-slip top copies interior u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(0, j, data->nx)], 0.0, 1e-10, "Free-slip bottom sets v=0");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(data->ny - 1, j, data->nx)], 0.0, 1e-10, "Free-slip top sets v=0");
    }

    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_NO_SLIP, 1.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_PERIODIC, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_RIGHT, BC_PERIODIC, 0.0, 0.0, NAN);

    memset(data->u, 0, (size_t)data->nx * data->ny * sizeof(double));
    memset(data->v, 0, (size_t)data->nx * data->ny * sizeof(double));
    memset(data->p, 0, (size_t)data->nx * data->ny * sizeof(double));

    for (int i = 0; i < data->ny; ++i) {
        data->u[IDX2D(i, 1, data->nx)] = 0.1 + 0.01 * i;
        data->v[IDX2D(i, 1, data->nx)] = 0.2;
        data->p[IDX2D(i, 1, data->nx)] = 1.0;

        data->u[IDX2D(i, data->nx - 2, data->nx)] = -0.2 + 0.02 * i;
        data->v[IDX2D(i, data->nx - 2, data->nx)] = -0.1;
        data->p[IDX2D(i, data->nx - 2, data->nx)] = 0.5;
    }

    ns_fd_apply_boundary_conditions(data);

    for (int i = 0; i < data->ny; ++i) {
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(i, 0, data->nx)], data->u[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Periodic left matches opposite interior u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(i, 0, data->nx)], data->v[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Periodic left matches opposite interior v");
        TEST_ASSERT_DOUBLE_EQUAL(data->p[IDX2D(i, 0, data->nx)], data->p[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Periodic left matches opposite interior p");

        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(i, data->nx - 1, data->nx)], data->u[IDX2D(i, 1, data->nx)], 1e-10, "Periodic right matches opposite interior u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(i, data->nx - 1, data->nx)], data->v[IDX2D(i, 1, data->nx)], 1e-10, "Periodic right matches opposite interior v");
        TEST_ASSERT_DOUBLE_EQUAL(data->p[IDX2D(i, data->nx - 1, data->nx)], data->p[IDX2D(i, 1, data->nx)], 1e-10, "Periodic right matches opposite interior p");
    }

    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_RIGHT, BC_OUTFLOW, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_INFLOW, 2.0, 0.1, 101.0);
    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_NO_SLIP, 0.0, 0.0, NAN);

    memset(data->u, 0, (size_t)data->nx * data->ny * sizeof(double));
    memset(data->v, 0, (size_t)data->nx * data->ny * sizeof(double));
    memset(data->p, 0, (size_t)data->nx * data->ny * sizeof(double));

    for (int i = 0; i < data->ny; ++i) {
        data->u[IDX2D(i, data->nx - 2, data->nx)] = 0.6 + 0.01 * i;
        data->v[IDX2D(i, data->nx - 2, data->nx)] = -0.15;
        data->p[IDX2D(i, data->nx - 2, data->nx)] = 4.0;
    }

    ns_fd_apply_boundary_conditions(data);

    for (int j = 0; j < data->nx; ++j) {
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(data->ny - 1, j, data->nx)], 2.0, 1e-10, "Inflow top enforces u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(data->ny - 1, j, data->nx)], 0.1, 1e-10, "Inflow top enforces v");
        TEST_ASSERT_DOUBLE_EQUAL(data->p[IDX2D(data->ny - 1, j, data->nx)], 101.0, 1e-10, "Inflow top enforces p");
    }

    for (int i = 0; i < data->ny; ++i) {
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(i, data->nx - 1, data->nx)], data->u[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Outflow right copies interior u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(i, data->nx - 1, data->nx)], data->v[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Outflow right copies interior v");
        TEST_ASSERT_DOUBLE_EQUAL(data->p[IDX2D(i, data->nx - 1, data->nx)], data->p[IDX2D(i, data->nx - 2, data->nx)], 1e-10, "Outflow right copies interior p");
    }

    CustomBCData custom = { -0.25, 0.5, 5.0 };
    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_CUSTOM, 0.0, 0.0, NAN);
    ns_fd_set_boundary_function(data, NS_BC_LEFT, custom_bc_callback, &custom);

    ns_fd_apply_boundary_conditions(data);

    for (int i = 0; i < data->ny; ++i) {
        TEST_ASSERT_DOUBLE_EQUAL(data->u[IDX2D(i, 0, data->nx)], custom.u, 1e-10, "Custom BC sets left u");
        TEST_ASSERT_DOUBLE_EQUAL(data->v[IDX2D(i, 0, data->nx)], custom.v, 1e-10, "Custom BC sets left v");
        TEST_ASSERT_DOUBLE_EQUAL(data->p[IDX2D(i, 0, data->nx)], custom.p, 1e-10, "Custom BC sets left p");
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

int test_common_forcing() {
    printf("  Testing turbulence forcing generator...\n");

    const int nkx = 4;
    const int nky = 4;
    double k2[nkx * nky];
    for (int j = 0; j < nky; ++j) {
        for (int i = 0; i < nkx; ++i) {
            const int idx = j * nkx + i;
            k2[idx] = (double)(i * i + j * j);
        }
    }

    NSForcingConfig cfg;
    ns_forcing_config_default(&cfg);
    cfg.kmin = 0.5;
    cfg.kmax = 2.5;
    cfg.energy_target = 1.5;
    cfg.seed = 42;
    cfg.type = NS_FORCING_WHITE_NOISE;

    NSForcingState state;
    int init_status = ns_forcing_state_init(&state, &cfg, k2, nkx, nky);
    TEST_ASSERT_EQUAL(init_status, 0, "Forcing state initialised");
    TEST_ASSERT(state.active_modes > 0, "Low-wavenumber modes selected");

    double complex forcing[nkx * nky];
    ns_forcing_generate(&state, 0.1, k2, forcing);

    double energy = 0.0;
    for (int idx = 0; idx < nkx * nky; ++idx) {
        double mag2 = creal(forcing[idx]) * creal(forcing[idx]) +
                      cimag(forcing[idx]) * cimag(forcing[idx]);
        energy += mag2;
    }
    TEST_ASSERT_DOUBLE_EQUAL(energy, cfg.energy_target, 1e-6, "Energy injection matches target");

    ns_forcing_generate(&state, 0.1, k2, forcing); // ensure repeatable
    ns_forcing_state_free(&state);

    cfg.type = NS_FORCING_COLORED_NOISE;
    cfg.correlation_time = 0.2;
    init_status = ns_forcing_state_init(&state, &cfg, k2, nkx, nky);
    TEST_ASSERT_EQUAL(init_status, 0, "Colored forcing state initialised");

    ns_forcing_generate(&state, 0.05, k2, forcing);
    double colored_energy = 0.0;
    for (int idx = 0; idx < nkx * nky; ++idx) {
        double mag2 = creal(forcing[idx]) * creal(forcing[idx]) +
                      cimag(forcing[idx]) * cimag(forcing[idx]);
        colored_energy += mag2;
    }
    TEST_ASSERT_DOUBLE_EQUAL(colored_energy, cfg.energy_target, 1e-6, "Colored forcing energy matches target");

    ns_forcing_state_free(&state);
    return 1;
}

int test_config_loader() {
    printf("  Testing configuration loader...\n");

    NSConfig config;
    ns_config_init(&config);

    TEST_ASSERT_EQUAL(ns_config_set(&config, "solver.type", "spectral"), 0, "Default solver set");
    TEST_ASSERT_EQUAL(ns_config_set(&config, "grid.nx", "64"), 0, "Default nx set");
    TEST_ASSERT_EQUAL(ns_config_set(&config, "grid.ny", "64"), 0, "Default ny set");
    TEST_ASSERT_EQUAL(ns_config_set(&config, "time.dt", "0.005"), 0, "Default dt set");

    const char *files[] = {"config/parameters.ini", "tests/config_override.json"};
    char error[256] = {0};
    int status = ns_config_load_files(&config, files, 2, error, sizeof(error));
    TEST_ASSERT_EQUAL(status, 0, error[0] ? error : "Configuration files loaded");

    TEST_ASSERT(ns_config_has_key(&config, "domain.Lx"), "Domain length present");
    TEST_ASSERT(strcmp(ns_config_get_string(&config, "solver.type", ""), "finite_difference") == 0,
                "Solver override applied");
    TEST_ASSERT_EQUAL(ns_config_get_int(&config, "grid.nx", 0), 192, "Grid nx override applied");
    TEST_ASSERT_EQUAL(ns_config_get_int(&config, "grid.ny", 0), 96, "Grid ny override applied");
    TEST_ASSERT_DOUBLE_EQUAL(ns_config_get_double(&config, "domain.Lx", 0.0), 2.0, 1e-12, "Domain Lx parsed");
    TEST_ASSERT_DOUBLE_EQUAL(ns_config_get_double(&config, "domain.Ly", 0.0), 2.0, 1e-12, "Domain Ly parsed");
    TEST_ASSERT_DOUBLE_EQUAL(ns_config_get_double(&config, "time.dt", 0.0), 0.001, 1e-12, "Time step override applied");
    TEST_ASSERT_DOUBLE_EQUAL(ns_config_get_double(&config, "forcing.energy_target", 0.0), 1.25, 1e-12,
                             "Forcing energy parsed");

    const char *solver_allowed[] = {"spectral", "finite_difference"};
    const char *forcing_allowed[] = {"none", "white", "colored"};

    NSConfigSchemaEntry schema[] = {
        {"solver.type", NS_CONFIG_VALUE_STRING, 1, 0.0, 0.0, 0, 0, solver_allowed, 2},
        {"grid.nx", NS_CONFIG_VALUE_INT, 1, 2.0, 0.0, 1, 0, NULL, 0},
        {"grid.ny", NS_CONFIG_VALUE_INT, 1, 2.0, 0.0, 1, 0, NULL, 0},
        {"domain.Lx", NS_CONFIG_VALUE_DOUBLE, 1, 1e-6, 0.0, 1, 0, NULL, 0},
        {"domain.Ly", NS_CONFIG_VALUE_DOUBLE, 1, 1e-6, 0.0, 1, 0, NULL, 0},
        {"fluid.Re", NS_CONFIG_VALUE_DOUBLE, 1, 1e-6, 0.0, 1, 0, NULL, 0},
        {"time.dt", NS_CONFIG_VALUE_DOUBLE, 1, 1e-9, 0.0, 1, 0, NULL, 0},
        {"time.t_end", NS_CONFIG_VALUE_DOUBLE, 1, 1e-9, 0.0, 1, 0, NULL, 0},
        {"forcing.energy_target", NS_CONFIG_VALUE_DOUBLE, 0, 0.0, 0.0, 1, 0, NULL, 0},
        {"forcing.type", NS_CONFIG_VALUE_STRING, 0, 0.0, 0.0, 0, 0, forcing_allowed, 3},
        {"spectral.hyperviscosity.coefficient", NS_CONFIG_VALUE_DOUBLE, 0, 0.0, 0.0, 1, 0, NULL, 0},
        {"spectral.hyperviscosity.order", NS_CONFIG_VALUE_INT, 0, 0.0, 0.0, 1, 0, NULL, 0}
    };

    status = ns_config_validate(&config, schema, sizeof(schema) / sizeof(schema[0]), error, sizeof(error));
    TEST_ASSERT_EQUAL(status, 0, error[0] ? error : "Configuration validated");

    ns_config_free(&config);
    return 1;
}

int test_fd_multigrid_poisson() {
    printf("  Testing FD multigrid Poisson solver...\n");

    const int nx = 33;
    const int ny = 33;
    const double L = 1.0;
    const double H = 1.0;
    const double pi = 3.14159265358979323846;

    NSFiniteDiffData *data = ns_fd_allocate(nx, ny, L, H, 100.0);
    TEST_ASSERT_NOT_NULL(data, "FD data allocated for multigrid test");

    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_DIRICHLET, 0.0, 0.0, 0.0);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_DIRICHLET, 0.0, 0.0, 0.0);
    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_DIRICHLET, 0.0, 0.0, 0.0);
    ns_fd_set_boundary_condition(data, NS_BC_RIGHT, BC_DIRICHLET, 0.0, 0.0, 0.0);

    const size_t total_points = (size_t)nx * (size_t)ny;
    double *rhs = (double *)malloc(total_points * sizeof(double));
    double *solution = (double *)calloc(total_points, sizeof(double));
    TEST_ASSERT_NOT_NULL(rhs, "RHS allocation successful");
    TEST_ASSERT_NOT_NULL(solution, "Solution allocation successful");

    for (int i = 0; i < ny; ++i) {
        const double y = i * data->dy;
        for (int j = 0; j < nx; ++j) {
            const double x = j * data->dx;
            const int idx = IDX2D(i, j, nx);
            if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1) {
                rhs[idx] = 0.0;
                solution[idx] = 0.0;
            } else {
                rhs[idx] = -2.0 * pi * pi * sin(pi * x) * sin(pi * y);
                solution[idx] = 0.0;
            }
        }
    }

    NSMultigridConfig config;
    ns_fd_multigrid_default_config(&config);
    config.tolerance = 1e-8;
    config.max_cycles = 60;
    config.pre_smoothing = 2;
    config.post_smoothing = 2;
    config.coarse_smoothing = 30;
    config.cycle_type = NS_MG_CYCLE_V;

    double history[128];
    NSMultigridMonitor monitor = { history, 128, 0 };

    NSMultigridStats stats = ns_fd_multigrid_solve_poisson(data, rhs, solution, &config, &monitor);

    TEST_ASSERT(stats.initial_residual > 0.0, "Initial residual computed");
    TEST_ASSERT(stats.converged, "Multigrid solver converged");
    TEST_ASSERT(stats.final_residual < 5e-7, "Final residual below tolerance");

    double max_error = 0.0;
    double l2_sum = 0.0;
    int interior_points = 0;

    for (int i = 1; i < ny - 1; ++i) {
        const double y = i * data->dy;
        for (int j = 1; j < nx - 1; ++j) {
            const double x = j * data->dx;
            const int idx = IDX2D(i, j, nx);
            const double exact = sin(pi * x) * sin(pi * y);
            const double diff = solution[idx] - exact;
            max_error = fmax(max_error, fabs(diff));
            l2_sum += diff * diff;
            ++interior_points;
        }
    }

    const double l2_error = sqrt(l2_sum / (double)interior_points);
    TEST_ASSERT(max_error < 5e-3, "Maximum error within tolerance");
    TEST_ASSERT(l2_error < 1e-3, "L2 error within tolerance");

    free(rhs);
    free(solution);
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

int test_spectral_hyperviscosity() {
    printf("  Testing spectral hyperviscosity configuration...\n");

    NSSpectralData *data = ns_spectral_allocate(32, 32, 2*M_PI, 2*M_PI, 100.0);
    TEST_ASSERT_NOT_NULL(data, "Spectral data allocated");

    ns_spectral_initialize_grids(data);

    NSSpectralHyperviscosityConfig cfg;
    ns_spectral_hyperviscosity_default(&cfg);
    cfg.enabled = 1;
    cfg.coefficient = 0.05;
    cfg.order = 4;
    cfg.cutoff_ratio = 0.6;
    cfg.sfd_enabled = 1;
    cfg.sfd_strength = 0.3;
    cfg.sfd_timescale = 0.4;

    ns_spectral_hyperviscosity_set(data, &cfg);

    double low_weight = data->hyperviscosity_weights[0];
    TEST_ASSERT_DOUBLE_EQUAL(low_weight, 0.0, 1e-12, "Low wavenumber weight zero");

    int high_idx = data->nkx * data->nky - 1;
    TEST_ASSERT(data->hyperviscosity_weights[high_idx] > 0.0, "High wavenumber weight positive");

    fftw_complex *rhs = data->rk1_hat;
    fftw_complex *state = data->omega_hat;

    for (int idx = 0; idx < data->nkx * data->nky; ++idx) {
        rhs[idx] = 0.0 + 0.0*I;
        state[idx] = 1.0 + 0.0*I;
    }

    ns_spectral_hyperviscosity_apply(data, rhs, state);

    double low_mag = hypot(creal(rhs[0]), cimag(rhs[0]));
    TEST_ASSERT(low_mag < 1e-8, "Low wavenumbers remain undamped");

    double high_mag = hypot(creal(rhs[high_idx]), cimag(rhs[high_idx]));
    TEST_ASSERT(high_mag > 1e-3, "High wavenumbers damped");

    ns_spectral_hyperviscosity_post_step(data, state, 0.1);
    double filtered = hypot(creal(data->sfd_hat[high_idx]), cimag(data->sfd_hat[high_idx]));
    TEST_ASSERT(filtered > 0.0, "Filtered state updated");

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
    RUN_TEST(test_fd_multigrid_poisson);
    RUN_TEST(test_common_forcing);
    RUN_TEST(test_config_loader);
    
    // Spectral tests
    RUN_TEST(test_spectral_allocation);
    RUN_TEST(test_spectral_grid_initialization);
    RUN_TEST(test_spectral_fft_transforms);
    RUN_TEST(test_spectral_dealiasing);
    RUN_TEST(test_spectral_velocity_from_vorticity);
    RUN_TEST(test_spectral_energy_computation);
    RUN_TEST(test_spectral_enstrophy_computation);
    RUN_TEST(test_spectral_hyperviscosity);
    
    TEST_SUITE_END();
}

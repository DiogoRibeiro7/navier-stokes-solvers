#ifndef NS_SPECTRAL_SOLVER_H
#define NS_SPECTRAL_SOLVER_H

#include "../../include/common_types.h"
#include <complex.h>
#include <fftw3.h>

// Spectral Navier-Stokes data structure
typedef struct {
    int nx, ny;                    // Grid dimensions
    int nkx, nky;                  // Fourier mode dimensions
    double Lx, Ly;                 // Domain size
    double Re;                     // Reynolds number
    double dt;                     // Time step
    double t;                      // Current time
    double cfl_factor;             // CFL safety factor
    
    // Physical space arrays
    double *u, *v, *p;            // Velocity and pressure
    double *omega;                // Vorticity
    double *psi;                  // Stream function
    double *x, *y;                // Grid coordinates
    
    // Spectral space arrays
    fftw_complex *u_hat, *v_hat;  // Velocity in Fourier space
    fftw_complex *omega_hat;      // Vorticity in Fourier space
    fftw_complex *psi_hat;        // Stream function
    fftw_complex *nonlinear_hat;  // Nonlinear terms
    
    // Wavenumber arrays
    double *kx, *ky;              // Wavenumber grids
    double *k2;                   // k^2 for Laplacian
    
    // FFT plans
    fftw_plan forward_plan, backward_plan;
    
    // Dealiasing mask
    int *dealias_mask;
    
} NSSpectralData;

// Memory management
NSSpectralData* ns_spectral_allocate(int nx, int ny, double Lx, double Ly, double Re);
void ns_spectral_free(NSSpectralData *data);

// Initialization
void ns_spectral_initialize_grids(NSSpectralData *data);
void ns_spectral_initialize_fft_plans(NSSpectralData *data);
void ns_spectral_initialize_taylor_green(NSSpectralData *data);

// Transforms
void ns_spectral_physical_to_spectral(NSSpectralData *data, double *field, fftw_complex *field_hat);
void ns_spectral_spectral_to_physical(NSSpectralData *data, fftw_complex *field_hat, double *field);
void ns_spectral_apply_dealiasing(NSSpectralData *data, fftw_complex *field_hat);

// Core solver functions
void ns_spectral_compute_derivatives(NSSpectralData *data);
void ns_spectral_compute_velocity_from_vorticity(NSSpectralData *data);
void ns_spectral_compute_nonlinear_terms(NSSpectralData *data);

// Time integration
void ns_spectral_advance_rk4(NSSpectralData *data);
double ns_spectral_compute_cfl_timestep(NSSpectralData *data);

// Analysis
double ns_spectral_compute_energy(NSSpectralData *data);
double ns_spectral_compute_enstrophy(NSSpectralData *data);
void ns_spectral_analyze_convergence(NSSpectralData *data, int step);
void ns_spectral_analyze_spectrum(NSSpectralData *data);
int ns_spectral_check_resolution(NSSpectralData *data);

// Output
void ns_spectral_output_solution(NSSpectralData *data, const char *filename);

#endif // NS_SPECTRAL_SOLVER_H

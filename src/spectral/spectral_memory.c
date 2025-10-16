#include "ns_spectral_solver.h"

// Allocate spectral solver data
NSSpectralData* ns_spectral_allocate(int nx, int ny, double Lx, double Ly, double Re) {
    NSSpectralData *data = malloc(sizeof(NSSpectralData));
    if (!data) return NULL;
    
    data->nx = nx;
    data->ny = ny;
    data->nkx = nx/2 + 1;
    data->nky = ny;
    data->Lx = Lx;
    data->Ly = Ly;
    data->Re = Re;
    data->dt = 0.01;
    data->t = 0.0;
    data->cfl_factor = 0.3;
    
    // Physical space
    data->u = fftw_alloc_real(nx * ny);
    data->v = fftw_alloc_real(nx * ny);
    data->p = fftw_alloc_real(nx * ny);
    data->omega = fftw_alloc_real(nx * ny);
    data->psi = fftw_alloc_real(nx * ny);
    data->x = malloc(nx * sizeof(double));
    data->y = malloc(ny * sizeof(double));
    
    // Spectral space
    int spectral_size = data->nkx * data->nky;
    data->u_hat = fftw_alloc_complex(spectral_size);
    data->v_hat = fftw_alloc_complex(spectral_size);
    data->omega_hat = fftw_alloc_complex(spectral_size);
    data->psi_hat = fftw_alloc_complex(spectral_size);
    data->nonlinear_hat = fftw_alloc_complex(spectral_size);
    
    // Wavenumbers
    data->kx = malloc(data->nkx * sizeof(double));
    data->ky = malloc(data->nky * sizeof(double));
    data->k2 = malloc(spectral_size * sizeof(double));
    data->dealias_mask = malloc(spectral_size * sizeof(int));
    
    // Check allocations
    if (!data->u || !data->v || !data->p || !data->omega || !data->psi ||
        !data->u_hat || !data->v_hat || !data->omega_hat || !data->psi_hat ||
        !data->nonlinear_hat || !data->kx || !data->ky || !data->k2 || 
        !data->dealias_mask || !data->x || !data->y) {
        ns_spectral_free(data);
        return NULL;
    }
    
    return data;
}

// Free spectral solver data
void ns_spectral_free(NSSpectralData *data) {
    if (!data) return;
    
    if (data->forward_plan) fftw_destroy_plan(data->forward_plan);
    if (data->backward_plan) fftw_destroy_plan(data->backward_plan);
    
    if (data->u) fftw_free(data->u);
    if (data->v) fftw_free(data->v);
    if (data->p) fftw_free(data->p);
    if (data->omega) fftw_free(data->omega);
    if (data->psi) fftw_free(data->psi);
    if (data->u_hat) fftw_free(data->u_hat);
    if (data->v_hat) fftw_free(data->v_hat);
    if (data->omega_hat) fftw_free(data->omega_hat);
    if (data->psi_hat) fftw_free(data->psi_hat);
    if (data->nonlinear_hat) fftw_free(data->nonlinear_hat);
    
    free(data->x); free(data->y);
    free(data->kx); free(data->ky); free(data->k2);
    free(data->dealias_mask);
    free(data);
    
    fftw_cleanup();
}

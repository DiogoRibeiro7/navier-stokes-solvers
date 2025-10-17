#include "ns_spectral_solver.h"
#include "../common/performance.h"
#include <string.h>

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
    data->forward_plan = NULL;
    data->backward_plan = NULL;
    
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
    data->rk1_hat = fftw_alloc_complex(spectral_size);
    data->rk2_hat = fftw_alloc_complex(spectral_size);
    data->rk3_hat = fftw_alloc_complex(spectral_size);
    data->rk4_hat = fftw_alloc_complex(spectral_size);
    data->omega_tmp_hat = fftw_alloc_complex(spectral_size);
    data->sfd_hat = fftw_alloc_complex(spectral_size);
    
    // Wavenumbers
    data->kx = malloc(data->nkx * sizeof(double));
    data->ky = malloc(data->nky * sizeof(double));
    data->k2 = malloc(spectral_size * sizeof(double));
    data->dealias_mask = malloc(spectral_size * sizeof(int));
    data->hyperviscosity_weights = malloc(spectral_size * sizeof(double));
    
    // Check allocations
    if (!data->u || !data->v || !data->p || !data->omega || !data->psi ||
        !data->u_hat || !data->v_hat || !data->omega_hat || !data->psi_hat ||
        !data->nonlinear_hat || !data->rk1_hat || !data->rk2_hat || !data->rk3_hat ||
        !data->rk4_hat || !data->omega_tmp_hat || !data->sfd_hat || !data->kx || !data->ky || !data->k2 || 
        !data->dealias_mask || !data->x || !data->y || !data->hyperviscosity_weights) {
        ns_spectral_free(data);
        return NULL;
    }

    ns_spectral_hyperviscosity_default(&data->hyperviscosity);
    memset(data->hyperviscosity_weights, 0, spectral_size * sizeof(double));
    memset(data->sfd_hat, 0, spectral_size * sizeof(fftw_complex));
    
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
    if (data->rk1_hat) fftw_free(data->rk1_hat);
    if (data->rk2_hat) fftw_free(data->rk2_hat);
    if (data->rk3_hat) fftw_free(data->rk3_hat);
    if (data->rk4_hat) fftw_free(data->rk4_hat);
    if (data->omega_tmp_hat) fftw_free(data->omega_tmp_hat);
    if (data->sfd_hat) fftw_free(data->sfd_hat);
    
    free(data->x); free(data->y);
    free(data->kx); free(data->ky); free(data->k2);
    free(data->dealias_mask);
    free(data->hyperviscosity_weights);
    free(data);
    
    fftw_cleanup();
#ifdef _OPENMP
    fftw_cleanup_threads();
#endif
}

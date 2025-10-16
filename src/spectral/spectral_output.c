#include "ns_spectral_solver.h"

// Output solution to file
void ns_spectral_output_solution(NSSpectralData *data, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return;
    
    int nx = data->nx, ny = data->ny;
    
    fprintf(fp, "# Spectral Navier-Stokes Solution\n");
    fprintf(fp, "# x y u v omega psi\n");
    
    // Update stream function in physical space
    ns_spectral_spectral_to_physical(data, data->psi_hat, data->psi);
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;
            fprintf(fp, "%.8f %.8f %.8e %.8e %.8e %.8e\n",
                   data->x[i], data->y[j], data->u[idx], data->v[idx],
                   data->omega[idx], data->psi[idx]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

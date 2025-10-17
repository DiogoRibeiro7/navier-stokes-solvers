#ifndef NS_VTK_WRITER_H
#define NS_VTK_WRITER_H

#include "../finite_difference/ns_fd_solver.h"
#include "../spectral/ns_spectral_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

int ns_vtk_write_fd_legacy(const char *filename,
                           const NSFiniteDiffData *state,
                           double time);

int ns_vtk_write_fd_xml(const char *filename,
                        const NSFiniteDiffData *state,
                        double time);

int ns_vtk_write_spectral_legacy(const char *filename,
                                 const NSSpectralData *state,
                                 double time);

int ns_vtk_write_spectral_xml(const char *filename,
                              const NSSpectralData *state,
                              double time);

int ns_vtk_update_pvd(const char *pvd_filename,
                      const char *dataset_file,
                      double time);

#ifdef __cplusplus
}
#endif

#endif /* NS_VTK_WRITER_H */

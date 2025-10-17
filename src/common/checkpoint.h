#ifndef NS_CHECKPOINT_H
#define NS_CHECKPOINT_H

#include "../finite_difference/ns_fd_solver.h"
#include "../spectral/ns_spectral_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

int ns_checkpoint_write_fd(const char *filename, const NSFiniteDiffData *state);
int ns_checkpoint_read_fd(const char *filename, NSFiniteDiffData *state);

int ns_checkpoint_write_spectral(const char *filename, const NSSpectralData *state);
int ns_checkpoint_read_spectral(const char *filename, NSSpectralData *state);

#ifdef __cplusplus
}
#endif

#endif /* NS_CHECKPOINT_H */

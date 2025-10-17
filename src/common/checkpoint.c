#include "checkpoint.h"

#include <hdf5.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NS_HDF5_COMP_LEVEL 4
#define NS_HDF5_CHUNK_MIN  8

static int write_scalar_attr(hid_t loc, const char *name, hid_t h5_type, const void *value) {
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) return -1;
    hid_t attr = H5Acreate2(loc, name, h5_type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Sclose(space);
        return -1;
    }
    herr_t status = H5Awrite(attr, h5_type, value);
    H5Aclose(attr);
    H5Sclose(space);
    return status < 0 ? -1 : 0;
}

static hid_t create_compressed_dataset(hid_t parent,
                                       const char *name,
                                       int rank,
                                       const hsize_t *dims) {
    hid_t space = H5Screate_simple(rank, dims, NULL);
    if (space < 0) {
        return -1;
    }
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    if (plist < 0) {
        H5Sclose(space);
        return -1;
    }
    hsize_t chunk[3];
    for (int i = 0; i < rank; ++i) {
        const hsize_t dim = dims[i];
        chunk[i] = dim < NS_HDF5_CHUNK_MIN ? dim : NS_HDF5_CHUNK_MIN;
    }
    if (H5Pset_chunk(plist, rank, chunk) < 0) {
        H5Pclose(plist);
        H5Sclose(space);
        return -1;
    }
#ifdef H5_HAVE_FILTER_DEFLATE
    H5Pset_deflate(plist, NS_HDF5_COMP_LEVEL);
#endif
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(space);
    return dset;
}

static int write_dataset_2d(hid_t parent,
                            const char *name,
                            const double *data,
                            hsize_t dim0,
                            hsize_t dim1) {
    hsize_t dims[2] = {dim0, dim1};
    hid_t dset = create_compressed_dataset(parent, name, 2, dims);
    if (dset < 0) {
        return -1;
    }
    herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    return status < 0 ? -1 : 0;
}

static int read_dataset_2d(hid_t parent,
                           const char *name,
                           double *data,
                           hsize_t dim0,
                           hsize_t dim1) {
    hid_t dset = H5Dopen2(parent, name, H5P_DEFAULT);
    if (dset < 0) {
        return -1;
    }
    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        H5Dclose(dset);
        return -1;
    }
    int rank = H5Sget_simple_extent_ndims(space);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (rank != 2 || dims[0] != dim0 || dims[1] != dim1) {
        H5Sclose(space);
        H5Dclose(dset);
        return -1;
    }
    herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Sclose(space);
    H5Dclose(dset);
    return status < 0 ? -1 : 0;
}

static int write_dataset_complex(hid_t parent,
                                 const char *name,
                                 const fftw_complex *data,
                                 hsize_t dim0,
                                 hsize_t dim1) {
    hsize_t dims[3] = {dim0, dim1, 2};
    hid_t dset = create_compressed_dataset(parent, name, 3, dims);
    if (dset < 0) {
        return -1;
    }
    size_t total = (size_t)dim0 * (size_t)dim1;
    double *buffer = (double *)malloc(total * 2 * sizeof(double));
    if (!buffer) {
        H5Dclose(dset);
        return -1;
    }
    for (size_t idx = 0; idx < total; ++idx) {
        buffer[2 * idx] = creal(data[idx]);
        buffer[2 * idx + 1] = cimag(data[idx]);
    }
    herr_t status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    free(buffer);
    H5Dclose(dset);
    return status < 0 ? -1 : 0;
}

static int read_dataset_complex(hid_t parent,
                                const char *name,
                                fftw_complex *data,
                                hsize_t dim0,
                                hsize_t dim1) {
    hid_t dset = H5Dopen2(parent, name, H5P_DEFAULT);
    if (dset < 0) {
        return -1;
    }
    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        H5Dclose(dset);
        return -1;
    }
    hsize_t dims[3];
    int rank = H5Sget_simple_extent_ndims(space);
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (rank != 3 || dims[0] != dim0 || dims[1] != dim1 || dims[2] != 2) {
        H5Sclose(space);
        H5Dclose(dset);
        return -1;
    }
    size_t total = (size_t)dim0 * (size_t)dim1;
    double *buffer = (double *)malloc(total * 2 * sizeof(double));
    if (!buffer) {
        H5Sclose(space);
        H5Dclose(dset);
        return -1;
    }
    herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    if (status >= 0) {
        for (size_t idx = 0; idx < total; ++idx) {
            const double real_part = buffer[2 * idx];
            const double imag_part = buffer[2 * idx + 1];
            data[idx] = real_part + imag_part * I;
        }
    }
    free(buffer);
    H5Sclose(space);
    H5Dclose(dset);
    return status < 0 ? -1 : 0;
}

int ns_checkpoint_write_fd(const char *filename, const NSFiniteDiffData *state) {
    if (!filename || !state) return -1;
    int status = 0;
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return -1;

    const char solver_type[] = "finite_difference";
    status |= write_string_attr(file, "solver", solver_type);
    status |= write_scalar_attr(file, "time", H5T_NATIVE_DOUBLE, &state->t);
    status |= write_scalar_attr(file, "dt", H5T_NATIVE_DOUBLE, &state->dt);
    status |= write_scalar_attr(file, "Re", H5T_NATIVE_DOUBLE, &state->Re);
    status |= write_scalar_attr(file, "nx", H5T_NATIVE_INT, &state->nx);
    status |= write_scalar_attr(file, "ny", H5T_NATIVE_INT, &state->ny);
    status |= write_scalar_attr(file, "L", H5T_NATIVE_DOUBLE, &state->L);
    status |= write_scalar_attr(file, "H", H5T_NATIVE_DOUBLE, &state->H);

    hid_t group = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(file);
        return -1;
    }

    const hsize_t ny = (hsize_t)state->ny;
    const hsize_t nx = (hsize_t)state->nx;
    status |= write_dataset_2d(group, "u", state->u, ny, nx);
    status |= write_dataset_2d(group, "v", state->v, ny, nx);
    status |= write_dataset_2d(group, "p", state->p, ny, nx);
    status |= write_dataset_2d(group, "u_old", state->u_old, ny, nx);
    status |= write_dataset_2d(group, "v_old", state->v_old, ny, nx);
    status |= write_dataset_2d(group, "p_old", state->p_old, ny, nx);

    H5Gclose(group);
    H5Fclose(file);
    return status;
}

int ns_checkpoint_read_fd(const char *filename, NSFiniteDiffData *state) {
    if (!filename || !state) return -1;
    int status = 0;
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return -1;

    int nx = 0, ny = 0;
    double L = 0.0, H = 0.0, Re = 0.0, time = 0.0, dt = 0.0;

    hid_t attr = H5Aopen(file, "nx", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_INT, &nx); H5Aclose(attr); }
    attr = H5Aopen(file, "ny", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_INT, &ny); H5Aclose(attr); }
    attr = H5Aopen(file, "L", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &L); H5Aclose(attr); }
    attr = H5Aopen(file, "H", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &H); H5Aclose(attr); }
    attr = H5Aopen(file, "Re", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &Re); H5Aclose(attr); }
    attr = H5Aopen(file, "time", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &time); H5Aclose(attr); }
    attr = H5Aopen(file, "dt", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &dt); H5Aclose(attr); }

    if (nx != state->nx || ny != state->ny) {
        H5Fclose(file);
        return -1;
    }

    hid_t group = H5Gopen2(file, "/fields", H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(file);
        return -1;
    }

    const hsize_t dim0 = (hsize_t)state->ny;
    const hsize_t dim1 = (hsize_t)state->nx;
    status |= read_dataset_2d(group, "u", state->u, dim0, dim1);
    status |= read_dataset_2d(group, "v", state->v, dim0, dim1);
    status |= read_dataset_2d(group, "p", state->p, dim0, dim1);
    status |= read_dataset_2d(group, "u_old", state->u_old, dim0, dim1);
    status |= read_dataset_2d(group, "v_old", state->v_old, dim0, dim1);
    status |= read_dataset_2d(group, "p_old", state->p_old, dim0, dim1);

    H5Gclose(group);
    H5Fclose(file);

    if (status == 0) {
        state->L = L;
        state->H = H;
        state->Re = Re;
        state->t = time;
        state->dt = dt;
    }
    return status;
}

int ns_checkpoint_write_spectral(const char *filename, const NSSpectralData *state) {
    if (!filename || !state) return -1;
    int status = 0;
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return -1;

    const char solver_type[] = "spectral";
    status |= write_string_attr(file, "solver", solver_type);
    status |= write_scalar_attr(file, "time", H5T_NATIVE_DOUBLE, &state->t);
    status |= write_scalar_attr(file, "dt", H5T_NATIVE_DOUBLE, &state->dt);
    status |= write_scalar_attr(file, "Re", H5T_NATIVE_DOUBLE, &state->Re);
    status |= write_scalar_attr(file, "nx", H5T_NATIVE_INT, &state->nx);
    status |= write_scalar_attr(file, "ny", H5T_NATIVE_INT, &state->ny);
    status |= write_scalar_attr(file, "Lx", H5T_NATIVE_DOUBLE, &state->Lx);
    status |= write_scalar_attr(file, "Ly", H5T_NATIVE_DOUBLE, &state->Ly);

    hid_t group = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(file);
        return -1;
    }

    const hsize_t ny = (hsize_t)state->ny;
    const hsize_t nx = (hsize_t)state->nx;
    status |= write_dataset_2d(group, "u", state->u, ny, nx);
    status |= write_dataset_2d(group, "v", state->v, ny, nx);
    status |= write_dataset_2d(group, "omega", state->omega, ny, nx);
    status |= write_dataset_2d(group, "psi", state->psi, ny, nx);

    status |= write_dataset_complex(group, "u_hat", state->u_hat, state->nky, state->nkx);
    status |= write_dataset_complex(group, "v_hat", state->v_hat, state->nky, state->nkx);
    status |= write_dataset_complex(group, "omega_hat", state->omega_hat, state->nky, state->nkx);

    H5Gclose(group);
    H5Fclose(file);
    return status;
}

int ns_checkpoint_read_spectral(const char *filename, NSSpectralData *state) {
    if (!filename || !state) return -1;
    int status = 0;
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return -1;

    int nx = 0, ny = 0;
    double Lx = 0.0, Ly = 0.0, Re = 0.0, time = 0.0, dt = 0.0;
    hid_t attr = H5Aopen(file, "nx", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_INT, &nx); H5Aclose(attr); }
    attr = H5Aopen(file, "ny", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_INT, &ny); H5Aclose(attr); }
    attr = H5Aopen(file, "Lx", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &Lx); H5Aclose(attr); }
    attr = H5Aopen(file, "Ly", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &Ly); H5Aclose(attr); }
    attr = H5Aopen(file, "Re", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &Re); H5Aclose(attr); }
    attr = H5Aopen(file, "time", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &time); H5Aclose(attr); }
    attr = H5Aopen(file, "dt", H5P_DEFAULT);
    if (attr >= 0) { H5Aread(attr, H5T_NATIVE_DOUBLE, &dt); H5Aclose(attr); }

    if (nx != state->nx || ny != state->ny) {
        H5Fclose(file);
        return -1;
    }

    hid_t group = H5Gopen2(file, "/fields", H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(file);
        return -1;
    }

    const hsize_t dim0 = (hsize_t)state->ny;
    const hsize_t dim1 = (hsize_t)state->nx;
    status |= read_dataset_2d(group, "u", state->u, dim0, dim1);
    status |= read_dataset_2d(group, "v", state->v, dim0, dim1);
    status |= read_dataset_2d(group, "omega", state->omega, dim0, dim1);
    status |= read_dataset_2d(group, "psi", state->psi, dim0, dim1);

    status |= read_dataset_complex(group, "u_hat", state->u_hat, state->nky, state->nkx);
    status |= read_dataset_complex(group, "v_hat", state->v_hat, state->nky, state->nkx);
    status |= read_dataset_complex(group, "omega_hat", state->omega_hat, state->nky, state->nkx);

    H5Gclose(group);
    H5Fclose(file);

    if (status == 0) {
        state->Lx = Lx;
        state->Ly = Ly;
        state->Re = Re;
        state->t = time;
        state->dt = dt;
    }
    return status;
}
static int write_string_attr(hid_t loc, const char *name, const char *value) {
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) return -1;
    if (H5Tset_size(type, strlen(value) + 1) < 0) {
        H5Tclose(type);
        return -1;
    }
    int status = write_scalar_attr(loc, name, type, value);
    H5Tclose(type);
    return status;
}

#include "vtk_writer.h"

#include "../../include/common_types.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int write_legacy_structured(const char *filename,
                                   const double *x_coords,
                                   const double *y_coords,
                                   const double *u,
                                   const double *v,
                                   const double *scalar,
                                   int nx,
                                   int ny,
                                   const char *vector_name,
                                   const char *scalar_name,
                                   double time) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        return -1;
    }

    const int points = nx * ny;
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Navier-Stokes Solution t=%.6f\n", time);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", nx, ny);
    fprintf(fp, "POINTS %d double\n", points);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int idx = j * nx + i;
            fprintf(fp, "%.9f %.9f 0.0\n",
                    x_coords[i],
                    y_coords[j]);
        }
    }

    fprintf(fp, "FIELD FieldData 1\n");
    fprintf(fp, "TIME 1 1 double\n%.9f\n", time);

    fprintf(fp, "POINT_DATA %d\n", points);
    if (u && v && vector_name) {
        fprintf(fp, "VECTORS %s double\n", vector_name);
        for (int idx = 0; idx < points; ++idx) {
            fprintf(fp, "%.9f %.9f 0.0\n", u[idx], v[idx]);
        }
    }
    if (scalar && scalar_name) {
        fprintf(fp, "SCALARS %s double 1\n", scalar_name);
        fprintf(fp, "LOOKUP_TABLE default\n");
        for (int idx = 0; idx < points; ++idx) {
            fprintf(fp, "%.9f\n", scalar[idx]);
        }
    }
    fclose(fp);
    return 0;
}

static int write_rectilinear_xml(const char *filename,
                                 const double *x_coords,
                                 const double *y_coords,
                                 const double *u,
                                 const double *v,
                                 const double *scalar,
                                 int nx,
                                 int ny,
                                 const char *vector_name,
                                 const char *scalar_name,
                                 double time) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        return -1;
    }

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 0\">\n", nx - 1, ny - 1);
    fprintf(fp, "    <FieldData>\n");
    fprintf(fp, "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">%.9f</DataArray>\n", time);
    fprintf(fp, "    </FieldData>\n");
    fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n", nx - 1, ny - 1);
    fprintf(fp, "      <PointData");
    if (scalar && scalar_name) {
        fprintf(fp, " Scalars=\"%s\"", scalar_name);
    }
    if (u && v && vector_name) {
        fprintf(fp, " Vectors=\"%s\"", vector_name);
    }
    fprintf(fp, ">\n");

    if (scalar && scalar_name) {
        fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", scalar_name);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                fprintf(fp, " %.9f", scalar[j * nx + i]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "        </DataArray>\n");
    }

    if (u && v && vector_name) {
        fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", vector_name);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int idx = j * nx + i;
                fprintf(fp, " %.9f %.9f 0.0", u[idx], v[idx]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "        </DataArray>\n");
    }

    fprintf(fp, "      </PointData>\n");
    fprintf(fp, "      <CellData/>\n");
    fprintf(fp, "      <Coordinates>\n");

    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"X\" format=\"ascii\">\n");
    for (int i = 0; i < nx; ++i) {
        fprintf(fp, " %.9f", x_coords[i]);
    }
    fprintf(fp, "\n        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Y\" format=\"ascii\">\n");
    for (int j = 0; j < ny; ++j) {
        fprintf(fp, " %.9f", y_coords[j]);
    }
    fprintf(fp, "\n        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Z\" format=\"ascii\"> 0.0 </DataArray>\n");
    fprintf(fp, "      </Coordinates>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </RectilinearGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
    return 0;
}

static void build_fd_coordinates(const NSFiniteDiffData *state,
                                 double *x_coords,
                                 double *y_coords) {
    for (int i = 0; i < state->nx; ++i) {
        x_coords[i] = i * state->dx;
    }
    for (int j = 0; j < state->ny; ++j) {
        y_coords[j] = j * state->dy;
    }
}

int ns_vtk_write_fd_legacy(const char *filename,
                           const NSFiniteDiffData *state,
                           double time) {
    if (!filename || !state) return -1;
    double *x_coords = (double *)malloc((size_t)state->nx * sizeof(double));
    double *y_coords = (double *)malloc((size_t)state->ny * sizeof(double));
    if (!x_coords || !y_coords) {
        free(x_coords);
        free(y_coords);
        return -1;
    }
    build_fd_coordinates(state, x_coords, y_coords);
    int result = write_legacy_structured(filename,
                                         x_coords,
                                         y_coords,
                                         state->u,
                                         state->v,
                                         state->p,
                                         state->nx,
                                         state->ny,
                                         "velocity",
                                         "pressure",
                                         time);
    free(x_coords);
    free(y_coords);
    return result;
}

int ns_vtk_write_fd_xml(const char *filename,
                        const NSFiniteDiffData *state,
                        double time) {
    if (!filename || !state) return -1;
    double *x_coords = (double *)malloc((size_t)state->nx * sizeof(double));
    double *y_coords = (double *)malloc((size_t)state->ny * sizeof(double));
    if (!x_coords || !y_coords) {
        free(x_coords);
        free(y_coords);
        return -1;
    }
    build_fd_coordinates(state, x_coords, y_coords);
    int result = write_rectilinear_xml(filename,
                                       x_coords,
                                       y_coords,
                                       state->u,
                                       state->v,
                                       state->p,
                                       state->nx,
                                       state->ny,
                                       "velocity",
                                       "pressure",
                                       time);
    free(x_coords);
    free(y_coords);
    return result;
}

int ns_vtk_write_spectral_legacy(const char *filename,
                                 const NSSpectralData *state,
                                 double time) {
    if (!filename || !state) return -1;
    return write_legacy_structured(filename,
                                   state->x,
                                   state->y,
                                   state->u,
                                   state->v,
                                   state->omega,
                                   state->nx,
                                   state->ny,
                                   "velocity",
                                   "vorticity",
                                   time);
}

int ns_vtk_write_spectral_xml(const char *filename,
                              const NSSpectralData *state,
                              double time) {
    if (!filename || !state) return -1;
    return write_rectilinear_xml(filename,
                                 state->x,
                                 state->y,
                                 state->u,
                                 state->v,
                                 state->omega,
                                 state->nx,
                                 state->ny,
                                 "velocity",
                                 "vorticity",
                                 time);
}

static char *read_entire_file(const char *filename, long *length) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return NULL;
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return NULL;
    }
    long size = ftell(fp);
    if (size < 0) {
        fclose(fp);
        return NULL;
    }
    rewind(fp);
    char *buffer = (char *)malloc((size_t)size + 1);
    if (!buffer) {
        fclose(fp);
        return NULL;
    }
    size_t read = fread(buffer, 1, (size_t)size, fp);
    fclose(fp);
    buffer[read] = '\0';
    if (length) {
        *length = (long)read;
    }
    return buffer;
}

int ns_vtk_update_pvd(const char *pvd_filename,
                      const char *dataset_file,
                      double time) {
    if (!pvd_filename || !dataset_file) return -1;
    FILE *test = fopen(pvd_filename, "r");
    if (!test) {
        FILE *fp = fopen(pvd_filename, "w");
        if (!fp) return -1;
        fprintf(fp,
                "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                "  <Collection>\n"
                "    <DataSet timestep=\"%.9f\" group=\"\" part=\"0\" file=\"%s\"/>\n"
                "  </Collection>\n"
                "</VTKFile>\n",
                time,
                dataset_file);
        fclose(fp);
        return 0;
    }
    fclose(test);

    long existing_len = 0;
    char *existing = read_entire_file(pvd_filename, &existing_len);
    if (!existing) return -1;

    const char *needle = "</Collection>";
    char *pos = strstr(existing, needle);
    if (!pos) {
        free(existing);
        return -1;
    }

    size_t prefix_len = (size_t)(pos - existing);
    const char *suffix_start = pos;

    FILE *fp = fopen(pvd_filename, "w");
    if (!fp) {
        free(existing);
        return -1;
    }
    fwrite(existing, 1, prefix_len, fp);
    fprintf(fp, "    <DataSet timestep=\"%.9f\" group=\"\" part=\"0\" file=\"%s\"/>\n", time, dataset_file);
    fputs(suffix_start, fp);
    fclose(fp);
    free(existing);
    return 0;
}

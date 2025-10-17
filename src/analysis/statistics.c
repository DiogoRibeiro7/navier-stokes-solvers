#include "statistics.h"

#include <math.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void ns_stats_compute_two_point_corr_x(const double *fluc,
                                              int nx, int ny,
                                              int max_sep,
                                              double variance,
                                              double *corr) {
    if (!corr || max_sep <= 0) {
        return;
    }

    if (variance <= 0.0) {
        for (int s = 0; s < max_sep; ++s) {
            corr[s] = (s == 0) ? 1.0 : 0.0;
        }
        return;
    }

    const double norm = 1.0 / ((double)(nx * ny) * variance);

#pragma omp parallel for
    for (int s = 0; s < max_sep; ++s) {
        double sum = 0.0;
        for (int i = 0; i < ny; ++i) {
            const int row = i * nx;
            for (int j = 0; j < nx; ++j) {
                const int idx = row + j;
                const int jp = row + ((j + s) % nx);
                sum += fluc[idx] * fluc[jp];
            }
        }
        corr[s] = sum * norm;
    }
}

static void ns_stats_compute_two_point_corr_y(const double *fluc,
                                              int nx, int ny,
                                              int max_sep,
                                              double variance,
                                              double *corr) {
    if (!corr || max_sep <= 0) {
        return;
    }

    if (variance <= 0.0) {
        for (int s = 0; s < max_sep; ++s) {
            corr[s] = (s == 0) ? 1.0 : 0.0;
        }
        return;
    }

    const double norm = 1.0 / ((double)(nx * ny) * variance);

#pragma omp parallel for
    for (int s = 0; s < max_sep; ++s) {
        double sum = 0.0;
        for (int i = 0; i < ny; ++i) {
            const int row = i * nx;
            const int ip = ((i + s) % ny) * nx;
            for (int j = 0; j < nx; ++j) {
                const int idx = row + j;
                const int neighbor = ip + j;
                sum += fluc[idx] * fluc[neighbor];
            }
        }
        corr[s] = sum * norm;
    }
}

static double ns_stats_integral_scale_from_corr(const double *corr,
                                                int len,
                                                double spacing) {
    if (!corr || len < 2 || spacing <= 0.0) {
        return 0.0;
    }

    double integral = 0.0;
    for (int s = 0; s < len - 1; ++s) {
        const double c0 = corr[s];
        const double c1 = corr[s + 1];

        if (c0 <= 0.0) {
            break;
        }
        if (c1 <= 0.0) {
            const double frac = c0 / (c0 - c1);
            integral += 0.5 * c0 * frac * spacing;
            break;
        }

        integral += 0.5 * (c0 + c1) * spacing;
    }
    return integral;
}

static double ns_stats_compute_taylor_microscale(const double *fluc,
                                                 int nx, int ny,
                                                 double dx,
                                                 double variance) {
    if (!fluc || nx < 3 || variance <= 0.0 || dx <= 0.0) {
        return 0.0;
    }

    double derivative_sum = 0.0;

#pragma omp parallel for reduction(+:derivative_sum)
    for (int i = 0; i < ny; ++i) {
        const int row = i * nx;
        for (int j = 0; j < nx; ++j) {
            const int idx = row + j;
            const int jp = row + ((j + 1) % nx);
            const int jm = row + ((j - 1 + nx) % nx);
            const double du_dx = (fluc[jp] - fluc[jm]) / (2.0 * dx);
            derivative_sum += du_dx * du_dx;
        }
    }

    const double avg_derivative_sq = derivative_sum / (double)(nx * ny);
    if (avg_derivative_sq <= 0.0) {
        return 0.0;
    }
    return sqrt(variance / avg_derivative_sq);
}

static double ns_stats_compute_kolmogorov_scale(const double *u,
                                                const double *v,
                                                int nx, int ny,
                                                double dx, double dy,
                                                double nu) {
    if (!u || !v || nx < 2 || ny < 2 || dx <= 0.0 || dy <= 0.0 || nu <= 0.0) {
        return 0.0;
    }

    double strain_sum = 0.0;

#pragma omp parallel for reduction(+:strain_sum)
    for (int i = 0; i < ny; ++i) {
        const int row = i * nx;
        const int ip = ((i + 1) % ny) * nx;
        const int im = ((i - 1 + ny) % ny) * nx;
        for (int j = 0; j < nx; ++j) {
            const int idx = row + j;
            const int jp = row + ((j + 1) % nx);
            const int jm = row + ((j - 1 + nx) % nx);

            const double du_dx = (u[jp] - u[jm]) / (2.0 * dx);
            const double du_dy = (u[ip + j] - u[im + j]) / (2.0 * dy);
            const double dv_dx = (v[jp] - v[jm]) / (2.0 * dx);
            const double dv_dy = (v[ip + j] - v[im + j]) / (2.0 * dy);

            const double Sxx = du_dx;
            const double Syy = dv_dy;
            const double Sxy = 0.5 * (du_dy + dv_dx);
            const double invariant = Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy;
            strain_sum += invariant;
        }
    }

    const double avg_strain = strain_sum / (double)(nx * ny);
    const double epsilon = 2.0 * nu * avg_strain;

    if (epsilon <= 0.0) {
        return 0.0;
    }

    const double ratio = (nu * nu * nu) / epsilon;
    return pow(ratio, 0.25);
}

int ns_stats_compute_velocity_statistics(const double *NS_RESTRICT u,
                                         const double *NS_RESTRICT v,
                                         int nx, int ny,
                                         double dx, double dy,
                                         double nu,
                                         int max_sep_x, int max_sep_y,
                                         NSTurbulenceStatistics *stats,
                                         double *corr_x, double *corr_y) {
    if (!u || !v || !stats || nx <= 0 || ny <= 0 || dx <= 0.0 || dy <= 0.0 ||
        max_sep_x < 0 || max_sep_y < 0) {
        return -1;
    }

    const int total_points = nx * ny;
    const double inv_total = 1.0 / (double)total_points;

    double sum_u = 0.0, sum_v = 0.0;
#pragma omp parallel for reduction(+:sum_u, sum_v)
    for (int i = 0; i < ny; ++i) {
        const int row = i * nx;
        for (int j = 0; j < nx; ++j) {
            const int idx = row + j;
            sum_u += u[idx];
            sum_v += v[idx];
        }
    }

    stats->u.mean = sum_u * inv_total;
    stats->v.mean = sum_v * inv_total;

    double *u_fluc = (double *)malloc((size_t)total_points * sizeof(double));
    if (!u_fluc) {
        return -1;
    }
    double *v_fluc = (double *)malloc((size_t)total_points * sizeof(double));
    if (!v_fluc) {
        free(u_fluc);
        return -1;
    }

    double sum2_u = 0.0, sum3_u = 0.0, sum4_u = 0.0;
    double sum2_v = 0.0, sum3_v = 0.0, sum4_v = 0.0;

#pragma omp parallel for reduction(+:sum2_u, sum3_u, sum4_u, sum2_v, sum3_v, sum4_v)
    for (int i = 0; i < ny; ++i) {
        const int row = i * nx;
        for (int j = 0; j < nx; ++j) {
            const int idx = row + j;
            const double du = u[idx] - stats->u.mean;
            const double dv = v[idx] - stats->v.mean;
            u_fluc[idx] = du;
            v_fluc[idx] = dv;
            const double du2 = du * du;
            const double dv2 = dv * dv;
            sum2_u += du2;
            sum2_v += dv2;
            sum3_u += du2 * du;
            sum3_v += dv2 * dv;
            sum4_u += du2 * du2;
            sum4_v += dv2 * dv2;
        }
    }

    const double var_u = sum2_u * inv_total;
    const double var_v = sum2_v * inv_total;

    stats->u.rms = (var_u > 0.0) ? sqrt(var_u) : 0.0;
    stats->v.rms = (var_v > 0.0) ? sqrt(var_v) : 0.0;

    if (stats->u.rms > 0.0) {
        const double inv_rms3 = 1.0 / (stats->u.rms * stats->u.rms * stats->u.rms);
        const double inv_rms4 = inv_rms3 / stats->u.rms;
        stats->u.skewness = (sum3_u * inv_total) * inv_rms3;
        stats->u.kurtosis = (sum4_u * inv_total) * inv_rms4;
    } else {
        stats->u.skewness = 0.0;
        stats->u.kurtosis = 0.0;
    }

    if (stats->v.rms > 0.0) {
        const double inv_rms3 = 1.0 / (stats->v.rms * stats->v.rms * stats->v.rms);
        const double inv_rms4 = inv_rms3 / stats->v.rms;
        stats->v.skewness = (sum3_v * inv_total) * inv_rms3;
        stats->v.kurtosis = (sum4_v * inv_total) * inv_rms4;
    } else {
        stats->v.skewness = 0.0;
        stats->v.kurtosis = 0.0;
    }

    const int sep_x = (max_sep_x > nx) ? nx : max_sep_x;
    const int sep_y = (max_sep_y > ny) ? ny : max_sep_y;

    double *local_corr_x = NULL;
    double *local_corr_y = NULL;
    double *corr_x_buf = corr_x;
    double *corr_y_buf = corr_y;

    if (sep_x > 0) {
        if (!corr_x_buf) {
            local_corr_x = (double *)malloc((size_t)sep_x * sizeof(double));
            if (!local_corr_x) {
                free(v_fluc);
                free(u_fluc);
                return -1;
            }
            corr_x_buf = local_corr_x;
        }
        ns_stats_compute_two_point_corr_x(u_fluc, nx, ny, sep_x, var_u, corr_x_buf);
    } else if (corr_x_buf && max_sep_x > 0) {
        for (int s = 0; s < max_sep_x; ++s) {
            corr_x_buf[s] = (s == 0) ? 1.0 : 0.0;
        }
    }

    if (sep_y > 0) {
        if (!corr_y_buf) {
            local_corr_y = (double *)malloc((size_t)sep_y * sizeof(double));
            if (!local_corr_y) {
                free(local_corr_x);
                free(v_fluc);
                free(u_fluc);
                return -1;
            }
            corr_y_buf = local_corr_y;
        }
        ns_stats_compute_two_point_corr_y(u_fluc, nx, ny, sep_y, var_u, corr_y_buf);
    } else if (corr_y_buf && max_sep_y > 0) {
        for (int s = 0; s < max_sep_y; ++s) {
            corr_y_buf[s] = (s == 0) ? 1.0 : 0.0;
        }
    }

    if (corr_x && corr_x_buf == corr_x && sep_x < max_sep_x) {
        for (int s = sep_x; s < max_sep_x; ++s) {
            corr_x[s] = 0.0;
        }
    }

    if (corr_y && corr_y_buf == corr_y && sep_y < max_sep_y) {
        for (int s = sep_y; s < max_sep_y; ++s) {
            corr_y[s] = 0.0;
        }
    }

    stats->integral_length_scale_x = (var_u > 0.0 && sep_x > 1 && corr_x_buf)
        ? ns_stats_integral_scale_from_corr(corr_x_buf, sep_x, dx)
        : 0.0;

    stats->integral_length_scale_y = (var_u > 0.0 && sep_y > 1 && corr_y_buf)
        ? ns_stats_integral_scale_from_corr(corr_y_buf, sep_y, dy)
        : 0.0;

    stats->taylor_microscale = ns_stats_compute_taylor_microscale(u_fluc, nx, ny, dx, var_u);
    stats->kolmogorov_scale = ns_stats_compute_kolmogorov_scale(u, v, nx, ny, dx, dy, nu);

    free(local_corr_x);
    free(local_corr_y);

    free(v_fluc);
    free(u_fluc);

    return 0;
}

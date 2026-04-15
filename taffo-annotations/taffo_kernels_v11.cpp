/*
 * taffo_kernels_v11.cpp — Precision-tuning kernels for the Thomas algorithm.
 *
 * Compile with:
 *   taffo -O3 -std=gnu++20 -fopenmp -march=native \
 *         -Xdta --maxtotalbits=<BITS> \
 *         -c taffo-annotations/taffo_kernels_v11.cpp \
 *         -o scratch/taffo_kernels_v11_<BITS>.o
 *
 * Actual value ranges for the default test problem
 * (D=10000, dt=0.01, dx=20, init_dens=1000, decay=0.1):
 *
 *   lambda = D * dt / dx²        = 10000 * 0.01 / 400  = 0.25
 *   b[0] = b[n-1]                = 1 + decay*dt/dims + lambda        ≈ 1.251
 *   b[i] (interior)              = 1 + decay*dt/dims + 2*lambda      ≈ 1.501
 *   b'[i] = 1 / (b[i] - c² b'[i-1])                   ∈ [0.50, 0.72]
 *   c[s]                         = lambda               = 0.25        (positive)
 *   e[i] = c * b'[i-1]          = 0.25 * b'            ∈ [0.13, 0.18] ⊆ [0, 0.25]
 *
 *   density physical range       : [0, 1000]
 *   forward-sweep peak (uniform) : d_peak ≈ d0 / (1 - e_max)
 *                                         ≈ 1000 / (1 - 0.18) ≈ 1220
 *   conservative annotated range : [0, 1300]  (rounds to 11 integer bits)
 *
 * Bit budget for density range [0, 1300]:
 *   ceil(log2(1300)) = 11 integer bits (2^11 = 2048 > 1300)
 *   --maxtotalbits=32 → 21 fractional bits, step ≈ 5e-7
 *   --maxtotalbits=16 →  5 fractional bits, step ≈ 0.03
 *   --maxtotalbits= 8 → overflow (8 < 11 needed); TAFFO may saturate or error
 *
 * Memory layout (sxyz, s outermost):
 *   dens[s*nx*ny*nz + x*ny*nz + y*nz + z]
 *   b[s*n + i],  e[s*n + i]   (length n per substrate)
 *   c[s]                      (one scalar per substrate)
 *
 * MAX_NS = 16 substrates, MAX_N = 512 per dimension.
 */

#include "taffo_kernels_v11.h"

static constexpr int MAX_NS = 16;
static constexpr int MAX_N  = 512;

/* Annotation shorthands */
#define A_DENS   __attribute__((annotate("scalar(range(0,1300) final)")))
#define A_BPRIME __attribute__((annotate("scalar(range(0.5,1.0) final)")))
#define A_EC     __attribute__((annotate("scalar(range(0,0.25) final)")))
#define A_OUT    __attribute__((annotate("scalar(disabled)")))
#define A_TARGET __attribute__((annotate("target('result') scalar(range(0,1300) final)")))

/* ================================================================ */
/* x-dimension — 3D/2D  (one yz slice, M = ny*nz or ny)            */
/* ================================================================ */

void taffo_x_slice_v11(double* __restrict__ dens,
                       const double* __restrict__ b,
                       const double* __restrict__ c,
                       const double* __restrict__ e,
                       int ns, int nx, int M, int yz)
{
	float local_dens[MAX_NS * MAX_N] A_DENS;

	/* copy in */
	for (int s = 0; s < ns; s++)
		for (int i = 0; i < nx; i++)
			local_dens[s * nx + i] = (float)dens[s * nx * M + i * M + yz];

	/* forward sweep: d[i] += e[i-1] * d[i-1] */
	for (int i = 1; i < nx; i++) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * nx + i];
			double d_prv A_DENS;   d_prv = local_dens[s * nx + (i - 1)];
			double e_val A_EC;     e_val = e[s * nx + (i - 1)];
			local_dens[s * nx + i] = (float)(d_i + e_val * d_prv);
		}
	}

	/* scale last element */
	for (int s = 0; s < ns; s++) {
		double d_n A_DENS;   d_n = local_dens[s * nx + (nx - 1)];
		double b_n A_BPRIME; b_n = b[s * nx + (nx - 1)];
		local_dens[s * nx + (nx - 1)] = (float)(d_n * b_n);
	}

	/* backward substitution: d[i] = (d[i] + c[s]*d[i+1]) * b[i] */
	for (int i = nx - 2; i >= 0; i--) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * nx + i];
			double d_nxt A_DENS;   d_nxt = local_dens[s * nx + (i + 1)];
			double b_val A_BPRIME; b_val = b[s * nx + i];
			double c_val A_EC;     c_val = c[s];
			double result A_TARGET;
			result = (d_i + c_val * d_nxt) * b_val;
			local_dens[s * nx + i] = (float)result;
		}
	}

	/* copy out: scalar(disabled) forces TAFFO to insert dequantization */
	for (int s = 0; s < ns; s++)
		for (int i = 0; i < nx; i++) {
			double out A_OUT;
			out = local_dens[s * nx + i];
			dens[s * nx * M + i * M + yz] = out;
		}
}

/* ================================================================ */
/* x-dimension — 1D  (one substrate s)                              */
/* ================================================================ */

void taffo_x_slice_1d_v11(double* __restrict__ dens,
                           const double* __restrict__ b,
                           const double* __restrict__ c,
                           const double* __restrict__ e,
                           int ns, int nx, int s)
{
	float local_dens[MAX_N] A_DENS;

	for (int i = 0; i < nx; i++)
		local_dens[i] = (float)dens[s * nx + i];

	for (int i = 1; i < nx; i++) {
		double d_i   A_DENS;   d_i   = local_dens[i];
		double d_prv A_DENS;   d_prv = local_dens[i - 1];
		double e_val A_EC;     e_val = e[s * nx + (i - 1)];
		local_dens[i] = (float)(d_i + e_val * d_prv);
	}
	{
		double d_n A_DENS;   d_n = local_dens[nx - 1];
		double b_n A_BPRIME; b_n = b[s * nx + (nx - 1)];
		local_dens[nx - 1] = (float)(d_n * b_n);
	}
	for (int i = nx - 2; i >= 0; i--) {
		double d_i   A_DENS;   d_i   = local_dens[i];
		double d_nxt A_DENS;   d_nxt = local_dens[i + 1];
		double b_val A_BPRIME; b_val = b[s * nx + i];
		double c_val A_EC;     c_val = c[s];
		double result A_TARGET;
		result = (d_i + c_val * d_nxt) * b_val;
		local_dens[i] = (float)result;
	}

	for (int i = 0; i < nx; i++) {
		double out A_OUT;
		out = local_dens[i];
		dens[s * nx + i] = out;
	}
}

/* ================================================================ */
/* y-dimension — 3D  (one (z, x) pair)                             */
/* ================================================================ */

void taffo_y_row_3d_v11(double* __restrict__ dens,
                         const double* __restrict__ b,
                         const double* __restrict__ c,
                         const double* __restrict__ e,
                         int ns, int nx, int ny, int nz, int z, int x)
{
	float local_dens[MAX_NS * MAX_N] A_DENS;

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < ny; i++)
			local_dens[s * ny + i] = (float)dens[s * nx * ny * nz + x * ny * nz + i * nz + z];

	for (int i = 1; i < ny; i++) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * ny + i];
			double d_prv A_DENS;   d_prv = local_dens[s * ny + (i - 1)];
			double e_val A_EC;     e_val = e[s * ny + (i - 1)];
			local_dens[s * ny + i] = (float)(d_i + e_val * d_prv);
		}
	}
	for (int s = 0; s < ns; s++) {
		double d_n A_DENS;   d_n = local_dens[s * ny + (ny - 1)];
		double b_n A_BPRIME; b_n = b[s * ny + (ny - 1)];
		local_dens[s * ny + (ny - 1)] = (float)(d_n * b_n);
	}
	for (int i = ny - 2; i >= 0; i--) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * ny + i];
			double d_nxt A_DENS;   d_nxt = local_dens[s * ny + (i + 1)];
			double b_val A_BPRIME; b_val = b[s * ny + i];
			double c_val A_EC;     c_val = c[s];
			double result A_TARGET;
			result = (d_i + c_val * d_nxt) * b_val;
			local_dens[s * ny + i] = (float)result;
		}
	}

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < ny; i++) {
			double out A_OUT;
			out = local_dens[s * ny + i];
			dens[s * nx * ny * nz + x * ny * nz + i * nz + z] = out;
		}
}

/* ================================================================ */
/* y-dimension — 2D  (one x column)                                */
/* ================================================================ */

void taffo_y_col_2d_v11(double* __restrict__ dens,
                         const double* __restrict__ b,
                         const double* __restrict__ c,
                         const double* __restrict__ e,
                         int ns, int nx, int ny, int x)
{
	float local_dens[MAX_NS * MAX_N] A_DENS;

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < ny; i++)
			local_dens[s * ny + i] = (float)dens[s * nx * ny + x * ny + i];

	for (int i = 1; i < ny; i++) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * ny + i];
			double d_prv A_DENS;   d_prv = local_dens[s * ny + (i - 1)];
			double e_val A_EC;     e_val = e[s * ny + (i - 1)];
			local_dens[s * ny + i] = (float)(d_i + e_val * d_prv);
		}
	}
	for (int s = 0; s < ns; s++) {
		double d_n A_DENS;   d_n = local_dens[s * ny + (ny - 1)];
		double b_n A_BPRIME; b_n = b[s * ny + (ny - 1)];
		local_dens[s * ny + (ny - 1)] = (float)(d_n * b_n);
	}
	for (int i = ny - 2; i >= 0; i--) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * ny + i];
			double d_nxt A_DENS;   d_nxt = local_dens[s * ny + (i + 1)];
			double b_val A_BPRIME; b_val = b[s * ny + i];
			double c_val A_EC;     c_val = c[s];
			double result A_TARGET;
			result = (d_i + c_val * d_nxt) * b_val;
			local_dens[s * ny + i] = (float)result;
		}
	}

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < ny; i++) {
			double out A_OUT;
			out = local_dens[s * ny + i];
			dens[s * nx * ny + x * ny + i] = out;
		}
}

/* ================================================================ */
/* z-dimension — 3D  (one (y, x) pair)                             */
/* ================================================================ */

void taffo_z_row_v11(double* __restrict__ dens,
                     const double* __restrict__ b,
                     const double* __restrict__ c,
                     const double* __restrict__ e,
                     int ns, int nx, int ny, int nz, int y, int x)
{
	float local_dens[MAX_NS * MAX_N] A_DENS;

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < nz; i++)
			local_dens[s * nz + i] = (float)dens[s * nx * ny * nz + x * ny * nz + y * nz + i];

	for (int i = 1; i < nz; i++) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * nz + i];
			double d_prv A_DENS;   d_prv = local_dens[s * nz + (i - 1)];
			double e_val A_EC;     e_val = e[s * nz + (i - 1)];
			local_dens[s * nz + i] = (float)(d_i + e_val * d_prv);
		}
	}
	for (int s = 0; s < ns; s++) {
		double d_n A_DENS;   d_n = local_dens[s * nz + (nz - 1)];
		double b_n A_BPRIME; b_n = b[s * nz + (nz - 1)];
		local_dens[s * nz + (nz - 1)] = (float)(d_n * b_n);
	}
	for (int i = nz - 2; i >= 0; i--) {
		for (int s = 0; s < ns; s++) {
			double d_i   A_DENS;   d_i   = local_dens[s * nz + i];
			double d_nxt A_DENS;   d_nxt = local_dens[s * nz + (i + 1)];
			double b_val A_BPRIME; b_val = b[s * nz + i];
			double c_val A_EC;     c_val = c[s];
			double result A_TARGET;
			result = (d_i + c_val * d_nxt) * b_val;
			local_dens[s * nz + i] = (float)result;
		}
	}

	for (int s = 0; s < ns; s++)
		for (int i = 0; i < nz; i++) {
			double out A_OUT;
			out = local_dens[s * nz + i];
			dens[s * nx * ny * nz + x * ny * nz + y * nz + i] = out;
		}
}

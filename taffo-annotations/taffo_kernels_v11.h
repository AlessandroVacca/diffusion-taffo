#pragma once

/*
 * taffo_kernels_v11.h — Precision-tuning kernels for the Thomas algorithm.
 *
 * Same local-array TAFFO pattern as v5 (the only approach that produces
 * correct TAFFO conversion), but with accurate value ranges from actual
 * problem parameters (D=10000, dt=0.01, dx=20, init=1000).
 *
 * Compile at multiple precisions (controls fractional bits available):
 *
 *   # 32-bit fixpt: 21 fractional bits, step ≈ 5e-7
 *   taffo -O3 -std=gnu++20 -fopenmp -march=native -Xdta --maxtotalbits=32 \
 *         -c taffo-annotations/taffo_kernels_v11.cpp \
 *         -o scratch/taffo_kernels_v11_32.o
 *
 *   # 16-bit fixpt: 5 fractional bits, step ≈ 0.03
 *   taffo -O3 -std=gnu++20 -fopenmp -march=native -Xdta --maxtotalbits=16 \
 *         -c taffo-annotations/taffo_kernels_v11.cpp \
 *         -o scratch/taffo_kernels_v11_16.o
 *
 *   # 8-bit fixpt: representation saturates (range needs 11 bits → 8 not enough)
 *   taffo -O3 -std=gnu++20 -fopenmp -march=native -Xdta --maxtotalbits=8 \
 *         -c taffo-annotations/taffo_kernels_v11.cpp \
 *         -o scratch/taffo_kernels_v11_8.o
 *
 * Actual value ranges (verified from precompute_values + Thomas sweep analysis):
 *   density, forward-sweep intermediates : [0, 1300]
 *     (peak ≈ d_init / (1 - e_max) ≈ 1000 / (1 - 0.24) ≈ 1316; always positive)
 *   b'[i]  (inv modified diagonal)       : [0.5, 1.0]
 *   e[i]   (= c * b'[i-1])              : [0, 0.25]   (positive; c > 0, b' > 0)
 *   c[s]   (off-diagonal coefficient)   : [0, 0.25]   (positive)
 *
 * Bit budget analysis (density range 0..1300, needs 11 integer bits since 2^11=2048):
 *   --maxtotalbits=32 → 32 - 11 = 21 fractional bits, step ≈ 5e-7
 *   --maxtotalbits=16 → 16 - 11 =  5 fractional bits, step ≈ 0.03
 *   --maxtotalbits= 8 → 8  - 11 = overflow; TAFFO will clamp or warn
 *
 * Function naming: all functions have suffix _v11 to distinguish from v5.
 * All functions use sxyz memory layout:
 *   dens[s*nx*ny*nz + x*ny*nz + y*nz + z]
 *   b[s*n + i],  e[s*n + i]   (n = nx/ny/nz)
 *   c[s]                      (scalar per substrate)
 */

/* x-dimension — 3D/2D: one yz-slice, M = ny*nz (3D) or ny (2D) */
void taffo_x_slice_v11(double* dens,
                       const double* b,
                       const double* c,
                       const double* e,
                       int ns, int nx, int M, int yz);

/* x-dimension — 1D: one substrate s */
void taffo_x_slice_1d_v11(double* dens,
                           const double* b,
                           const double* c,
                           const double* e,
                           int ns, int nx, int s);

/* y-dimension — 3D: one (z, x) pair */
void taffo_y_row_3d_v11(double* dens,
                         const double* b,
                         const double* c,
                         const double* e,
                         int ns, int nx, int ny, int nz, int z, int x);

/* y-dimension — 2D: one x column */
void taffo_y_col_2d_v11(double* dens,
                         const double* b,
                         const double* c,
                         const double* e,
                         int ns, int nx, int ny, int x);

/* z-dimension — 3D: one (y, x) pair */
void taffo_z_row_v11(double* dens,
                     const double* b,
                     const double* c,
                     const double* e,
                     int ns, int nx, int ny, int nz, int y, int x);

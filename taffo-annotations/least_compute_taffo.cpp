/*
 * least_compute_taffo.cpp
 *
 * Drop-in replacement for src/least_compute_thomas_solver.cpp.
 * Compiled with regular clang-18 (not taffo).
 *
 * For double: delegates hot Thomas loops to taffo_kernels_v11.o (TAFFO-compiled).
 * For float:  uses plain flat-index Thomas (generic_* helpers below).
 *
 * Memory layouts:
 *   dens[s*nx*ny*nz + x*ny*nz + y*nz + z]   (sxyz, s outermost)
 *   b[s*n + i],  e[s*n + (i-1)]             (s outermost, length n per substrate)
 *   c[s]                                     (one scalar per substrate)
 */

#include "least_compute_thomas_solver.h"
#include "perf_utils.h"
#include "taffo_kernels_v11.h"

#include <memory>
#include <noarr/structures_extended.hpp>

/* ====================================================================== */
/* Generic flat Thomas helpers — used by the float path                    */
/* ====================================================================== */

template <typename real_t>
static inline void generic_x_slice(real_t* __restrict__ dens,
                                    const real_t* __restrict__ b,
                                    const real_t* __restrict__ c,
                                    const real_t* __restrict__ e,
                                    int ns, int nx, int M, int yz)
{
	for (int i = 1; i < nx; i++)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * M + i       * M + yz];
			const real_t d_prev = dens[s * nx * M + (i - 1) * M + yz];
			const real_t e_val  = e[s * nx + (i - 1)];
			dens[s * nx * M + i * M + yz] = d_i + e_val * d_prev;
		}
	for (int s = 0; s < ns; s++)
		dens[s * nx * M + (nx - 1) * M + yz] *= b[s * nx + (nx - 1)];
	for (int i = nx - 2; i >= 0; i--)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * M + i       * M + yz];
			const real_t d_next = dens[s * nx * M + (i + 1) * M + yz];
			const real_t b_val  = b[s * nx + i];
			const real_t c_val  = c[s];
			dens[s * nx * M + i * M + yz] = (d_i + c_val * d_next) * b_val;
		}
}

template <typename real_t>
static inline void generic_x_slice_1d(real_t* __restrict__ dens,
                                       const real_t* __restrict__ b,
                                       const real_t* __restrict__ c,
                                       const real_t* __restrict__ e,
                                       int ns, int nx, int s)
{
	for (int i = 1; i < nx; i++) {
		const real_t d_i    = dens[s * nx + i];
		const real_t d_prev = dens[s * nx + (i - 1)];
		const real_t e_val  = e[s * nx + (i - 1)];
		dens[s * nx + i] = d_i + e_val * d_prev;
	}
	dens[s * nx + (nx - 1)] *= b[s * nx + (nx - 1)];
	for (int i = nx - 2; i >= 0; i--) {
		const real_t d_i    = dens[s * nx + i];
		const real_t d_next = dens[s * nx + (i + 1)];
		const real_t b_val  = b[s * nx + i];
		const real_t c_val  = c[s];
		dens[s * nx + i] = (d_i + c_val * d_next) * b_val;
	}
}

template <typename real_t>
static inline void generic_y_slice_3d(real_t* __restrict__ dens,
                                       const real_t* __restrict__ b,
                                       const real_t* __restrict__ c,
                                       const real_t* __restrict__ e,
                                       int ns, int nx, int ny, int nz, int z)
{
	for (int i = 1; i < ny; i++)
		for (int x = 0; x < nx; x++)
			for (int s = 0; s < ns; s++) {
				const real_t d_i    = dens[s * nx * ny * nz + x * ny * nz + i       * nz + z];
				const real_t d_prev = dens[s * nx * ny * nz + x * ny * nz + (i - 1) * nz + z];
				const real_t e_val  = e[s * ny + (i - 1)];
				dens[s * nx * ny * nz + x * ny * nz + i * nz + z] = d_i + e_val * d_prev;
			}
	for (int x = 0; x < nx; x++)
		for (int s = 0; s < ns; s++)
			dens[s * nx * ny * nz + x * ny * nz + (ny - 1) * nz + z] *= b[s * ny + (ny - 1)];
	for (int i = ny - 2; i >= 0; i--)
		for (int x = 0; x < nx; x++)
			for (int s = 0; s < ns; s++) {
				const real_t d_i    = dens[s * nx * ny * nz + x * ny * nz + i       * nz + z];
				const real_t d_next = dens[s * nx * ny * nz + x * ny * nz + (i + 1) * nz + z];
				const real_t b_val  = b[s * ny + i];
				const real_t c_val  = c[s];
				dens[s * nx * ny * nz + x * ny * nz + i * nz + z] = (d_i + c_val * d_next) * b_val;
			}
}

template <typename real_t>
static inline void generic_y_slice_2d(real_t* __restrict__ dens,
                                       const real_t* __restrict__ b,
                                       const real_t* __restrict__ c,
                                       const real_t* __restrict__ e,
                                       int ns, int nx, int ny, int x)
{
	for (int i = 1; i < ny; i++)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * ny + x * ny + i];
			const real_t d_prev = dens[s * nx * ny + x * ny + (i - 1)];
			const real_t e_val  = e[s * ny + (i - 1)];
			dens[s * nx * ny + x * ny + i] = d_i + e_val * d_prev;
		}
	for (int s = 0; s < ns; s++)
		dens[s * nx * ny + x * ny + (ny - 1)] *= b[s * ny + (ny - 1)];
	for (int i = ny - 2; i >= 0; i--)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * ny + x * ny + i];
			const real_t d_next = dens[s * nx * ny + x * ny + (i + 1)];
			const real_t b_val  = b[s * ny + i];
			const real_t c_val  = c[s];
			dens[s * nx * ny + x * ny + i] = (d_i + c_val * d_next) * b_val;
		}
}

template <typename real_t>
static inline void generic_z_slice(real_t* __restrict__ dens,
                                    const real_t* __restrict__ b,
                                    const real_t* __restrict__ c,
                                    const real_t* __restrict__ e,
                                    int ns, int nx, int ny, int nz, int y, int x)
{
	for (int i = 1; i < nz; i++)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * ny * nz + x * ny * nz + y * nz + i];
			const real_t d_prev = dens[s * nx * ny * nz + x * ny * nz + y * nz + (i - 1)];
			const real_t e_val  = e[s * nz + (i - 1)];
			dens[s * nx * ny * nz + x * ny * nz + y * nz + i] = d_i + e_val * d_prev;
		}
	for (int s = 0; s < ns; s++)
		dens[s * nx * ny * nz + x * ny * nz + y * nz + (nz - 1)] *= b[s * nz + (nz - 1)];
	for (int i = nz - 2; i >= 0; i--)
		for (int s = 0; s < ns; s++) {
			const real_t d_i    = dens[s * nx * ny * nz + x * ny * nz + y * nz + i];
			const real_t d_next = dens[s * nx * ny * nz + x * ny * nz + y * nz + (i + 1)];
			const real_t b_val  = b[s * nz + i];
			const real_t c_val  = c[s];
			dens[s * nx * ny * nz + x * ny * nz + y * nz + i] = (d_i + c_val * d_next) * b_val;
		}
}

/* ====================================================================== */
/* precompute_values — identical to original                               */
/* ====================================================================== */

template <typename real_t>
void least_compute_thomas_solver<real_t>::precompute_values(std::unique_ptr<real_t[]>& b,
                                                             std::unique_ptr<real_t[]>& c,
                                                             std::unique_ptr<real_t[]>& e,
                                                             index_t shape, index_t dims,
                                                             index_t n, index_t copies)
{
	b = std::make_unique<real_t[]>(n * this->problem_.substrates_count * copies);
	e = std::make_unique<real_t[]>((n - 1) * this->problem_.substrates_count * copies);
	c = std::make_unique<real_t[]>(this->problem_.substrates_count * copies);

	auto layout = noarr::scalar<real_t>() ^ noarr::vector<'s'>(this->problem_.substrates_count)
				  ^ noarr::vector<'x'>(copies) ^ noarr::vector<'i'>(n);

	auto b_diag = noarr::make_bag(layout, b.get());
	auto e_diag = noarr::make_bag(layout, e.get());

	for (index_t x = 0; x < copies; x++)
		for (index_t s = 0; s < this->problem_.substrates_count; s++)
			c[x * this->problem_.substrates_count + s] =
				-1 * -this->problem_.dt * this->problem_.diffusion_coefficients[s] / (shape * shape);

	{
		std::array<index_t, 2> indices = { 0, n - 1 };

		for (index_t i : indices)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < this->problem_.substrates_count; s++)
					b_diag.template at<'i', 'x', 's'>(i, x, s) =
						1 + this->problem_.decay_rates[s] * this->problem_.dt / dims
						+ this->problem_.dt * this->problem_.diffusion_coefficients[s] / (shape * shape);

		for (index_t i = 1; i < n - 1; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < this->problem_.substrates_count; s++)
					b_diag.template at<'i', 'x', 's'>(i, x, s) =
						1 + this->problem_.decay_rates[s] * this->problem_.dt / dims
						+ 2 * this->problem_.dt * this->problem_.diffusion_coefficients[s] / (shape * shape);
	}

	{
		for (index_t x = 0; x < copies; x++)
			for (index_t s = 0; s < this->problem_.substrates_count; s++)
				b_diag.template at<'i', 'x', 's'>(0, x, s) = 1 / b_diag.template at<'i', 'x', 's'>(0, x, s);

		for (index_t i = 1; i < n; i++)
			for (index_t x = 0; x < copies; x++)
				for (index_t s = 0; s < this->problem_.substrates_count; s++)
				{
					b_diag.template at<'i', 'x', 's'>(i, x, s) =
						1
						/ (b_diag.template at<'i', 'x', 's'>(i, x, s)
						   - c[x * this->problem_.substrates_count + s]
								 * c[x * this->problem_.substrates_count + s]
								 * b_diag.template at<'i', 'x', 's'>(i - 1, x, s));

					e_diag.template at<'i', 'x', 's'>(i - 1, x, s) =
						c[x * this->problem_.substrates_count + s]
						* b_diag.template at<'i', 'x', 's'>(i - 1, x, s);
				}
	}
}

template <typename real_t>
void least_compute_thomas_solver<real_t>::initialize()
{
	if (this->problem_.dims >= 1)
		precompute_values(bx_, cx_, ex_, this->problem_.dx, this->problem_.dims, this->problem_.nx, 1);
	if (this->problem_.dims >= 2)
		precompute_values(by_, cy_, ey_, this->problem_.dy, this->problem_.dims, this->problem_.ny, 1);
	if (this->problem_.dims >= 3)
		precompute_values(bz_, cz_, ez_, this->problem_.dz, this->problem_.dims, this->problem_.nz, 1);
}

template <typename real_t>
auto least_compute_thomas_solver<real_t>::get_diagonal_layout(const problem_t<index_t, real_t>& problem, index_t n)
{
	return noarr::scalar<real_t>() ^ noarr::vectors<'s', 'i'>(problem.substrates_count, n);
}

/* ====================================================================== */
/* solve_x                                                                 */
/* ====================================================================== */

template <typename real_t>
void least_compute_thomas_solver<real_t>::solve_x()
{
	const int ns = (int)this->problem_.substrates_count;
	const int nx = (int)this->problem_.nx;
	real_t* dens  = this->substrates_;
	const real_t* bx = bx_.get();
	const real_t* cx = cx_.get();
	const real_t* ex = ex_.get();

	if (this->problem_.dims == 1)
	{
		if constexpr (std::is_same_v<real_t, double>) {
#pragma omp parallel for schedule(static)
			for (int s = 0; s < ns; s++)
				taffo_x_slice_1d_v11(dens, bx, cx, ex, ns, nx, s);
		} else {
#pragma omp parallel for schedule(static)
			for (int s = 0; s < ns; s++)
				generic_x_slice_1d(dens, bx, cx, ex, ns, nx, s);
		}
	}
	else
	{
		const int M = (this->problem_.dims == 2) ? (int)this->problem_.ny
		                                          : (int)(this->problem_.ny * this->problem_.nz);
		if constexpr (std::is_same_v<real_t, double>) {
#pragma omp parallel for schedule(static)
			for (int yz = 0; yz < M; yz++)
				taffo_x_slice_v11(dens, bx, cx, ex, ns, nx, M, yz);
		} else {
#pragma omp parallel for schedule(static)
			for (int yz = 0; yz < M; yz++)
				generic_x_slice(dens, bx, cx, ex, ns, nx, M, yz);
		}
	}
}

/* ====================================================================== */
/* solve_y                                                                 */
/* ====================================================================== */

template <typename real_t>
void least_compute_thomas_solver<real_t>::solve_y()
{
	const int ns = (int)this->problem_.substrates_count;
	const int nx = (int)this->problem_.nx;
	const int ny = (int)this->problem_.ny;
	real_t* dens  = this->substrates_;
	const real_t* by = by_.get();
	const real_t* cy = cy_.get();
	const real_t* ey = ey_.get();

	if (this->problem_.dims == 2)
	{
		if constexpr (std::is_same_v<real_t, double>) {
#pragma omp parallel for schedule(static)
			for (int x = 0; x < nx; x++)
				taffo_y_col_2d_v11(dens, by, cy, ey, ns, nx, ny, x);
		} else {
#pragma omp parallel for schedule(static)
			for (int x = 0; x < nx; x++)
				generic_y_slice_2d(dens, by, cy, ey, ns, nx, ny, x);
		}
	}
	else if (this->problem_.dims == 3)
	{
		const int nz = (int)this->problem_.nz;
		if constexpr (std::is_same_v<real_t, double>) {
#pragma omp parallel for schedule(static) collapse(2)
			for (int z = 0; z < nz; z++)
				for (int x = 0; x < nx; x++)
					taffo_y_row_3d_v11(dens, by, cy, ey, ns, nx, ny, nz, z, x);
		} else {
#pragma omp parallel for schedule(static)
			for (int z = 0; z < nz; z++)
				generic_y_slice_3d(dens, by, cy, ey, ns, nx, ny, nz, z);
		}
	}
}

/* ====================================================================== */
/* solve_z                                                                 */
/* ====================================================================== */

template <typename real_t>
void least_compute_thomas_solver<real_t>::solve_z()
{
	const int ns = (int)this->problem_.substrates_count;
	const int nx = (int)this->problem_.nx;
	const int ny = (int)this->problem_.ny;
	const int nz = (int)this->problem_.nz;
	real_t* dens  = this->substrates_;
	const real_t* bz = bz_.get();
	const real_t* cz = cz_.get();
	const real_t* ez = ez_.get();

	if constexpr (std::is_same_v<real_t, double>) {
#pragma omp parallel for schedule(static) collapse(2)
		for (int y = 0; y < ny; y++)
			for (int x = 0; x < nx; x++)
				taffo_z_row_v11(dens, bz, cz, ez, ns, nx, ny, nz, y, x);
	} else {
#pragma omp parallel for schedule(static) collapse(2)
		for (int y = 0; y < ny; y++)
			for (int x = 0; x < nx; x++)
				generic_z_slice(dens, bz, cz, ez, ns, nx, ny, nz, y, x);
	}
}

/* ====================================================================== */
/* solve — full time-stepping loop                                         */
/* ====================================================================== */

template <typename real_t>
void least_compute_thomas_solver<real_t>::solve()
{
	const int ns = (int)this->problem_.substrates_count;
	const int nx = (int)this->problem_.nx;
	real_t* dens  = this->substrates_;
	const real_t* bx = bx_.get();
	const real_t* cx = cx_.get();
	const real_t* ex = ex_.get();

	if (this->problem_.dims == 1)
	{
#pragma omp parallel
		{
			perf_counter counter("lstc");
			for (index_t iter = 0; iter < this->problem_.iterations; iter++) {
				if constexpr (std::is_same_v<real_t, double>) {
#pragma omp for schedule(static) nowait
					for (int s = 0; s < ns; s++)
						taffo_x_slice_1d_v11(dens, bx, cx, ex, ns, nx, s);
				} else {
#pragma omp for schedule(static) nowait
					for (int s = 0; s < ns; s++)
						generic_x_slice_1d(dens, bx, cx, ex, ns, nx, s);
				}
#pragma omp barrier
			}
		}
	}
	else if (this->problem_.dims == 2)
	{
		const int ny = (int)this->problem_.ny;
		const int M  = ny;
		const real_t* by = by_.get();
		const real_t* cy = cy_.get();
		const real_t* ey = ey_.get();

#pragma omp parallel
		{
			perf_counter counter("lstc");
			for (index_t iter = 0; iter < this->problem_.iterations; iter++) {
				if constexpr (std::is_same_v<real_t, double>) {
#pragma omp for schedule(static) nowait
					for (int yz = 0; yz < M; yz++)
						taffo_x_slice_v11(dens, bx, cx, ex, ns, nx, M, yz);
#pragma omp barrier
#pragma omp for schedule(static) nowait
					for (int x = 0; x < nx; x++)
						taffo_y_col_2d_v11(dens, by, cy, ey, ns, nx, ny, x);
#pragma omp barrier
				} else {
#pragma omp for schedule(static) nowait
					for (int yz = 0; yz < M; yz++)
						generic_x_slice(dens, bx, cx, ex, ns, nx, M, yz);
#pragma omp barrier
#pragma omp for schedule(static) nowait
					for (int x = 0; x < nx; x++)
						generic_y_slice_2d(dens, by, cy, ey, ns, nx, ny, x);
#pragma omp barrier
				}
			}
		}
	}
	else if (this->problem_.dims == 3)
	{
		const int ny = (int)this->problem_.ny;
		const int nz = (int)this->problem_.nz;
		const int M  = ny * nz;
		const real_t* by = by_.get();
		const real_t* cy = cy_.get();
		const real_t* ey = ey_.get();
		const real_t* bz = bz_.get();
		const real_t* cz = cz_.get();
		const real_t* ez = ez_.get();

#pragma omp parallel
		{
			perf_counter counter("lstc");
			for (index_t iter = 0; iter < this->problem_.iterations; iter++) {
				if constexpr (std::is_same_v<real_t, double>) {
#pragma omp for schedule(static) nowait
					for (int yz = 0; yz < M; yz++)
						taffo_x_slice_v11(dens, bx, cx, ex, ns, nx, M, yz);
#pragma omp barrier
#pragma omp for schedule(static) collapse(2) nowait
					for (int z = 0; z < nz; z++)
						for (int x = 0; x < nx; x++)
							taffo_y_row_3d_v11(dens, by, cy, ey, ns, nx, ny, nz, z, x);
#pragma omp barrier
#pragma omp for schedule(static) collapse(2) nowait
					for (int y = 0; y < ny; y++)
						for (int x = 0; x < nx; x++)
							taffo_z_row_v11(dens, bz, cz, ez, ns, nx, ny, nz, y, x);
#pragma omp barrier
				} else {
#pragma omp for schedule(static) nowait
					for (int yz = 0; yz < M; yz++)
						generic_x_slice(dens, bx, cx, ex, ns, nx, M, yz);
#pragma omp barrier
#pragma omp for schedule(static) nowait
					for (int z = 0; z < nz; z++)
						generic_y_slice_3d(dens, by, cy, ey, ns, nx, ny, nz, z);
#pragma omp barrier
#pragma omp for schedule(static) collapse(2) nowait
					for (int y = 0; y < ny; y++)
						for (int x = 0; x < nx; x++)
							generic_z_slice(dens, bz, cz, ez, ns, nx, ny, nz, y, x);
#pragma omp barrier
				}
			}
		}
	}
}

/* ====================================================================== */
/* Explicit instantiations                                                 */
/* ====================================================================== */

template class least_compute_thomas_solver<float>;
template class least_compute_thomas_solver<double>;

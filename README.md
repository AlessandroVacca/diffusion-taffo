# Port of diffusion-benchmarking to fixed point math
This repo is the port of the `lstc` (Least Compute Thomas) diffusion solver in this [repository](https://github.com/DARE-PhysiBoSS/diffusion-benchmarking) to fixed point math, by using [TAFFO](https://github.com/TAFFO-org/TAFFO). The goal was to measure
the performance and accuracy tradeoff of replacing double-precision arithmetic with reduced-precision
fixed-point arithmetic in the various numerical kernel of the solver.

The baseline for all comparisons is the original `lstc` solver compiled with GCC,
`-march=native`, and `--double` (double precision). The annotated version (`diffuse-taffo`)
is compiled with clang-18 and always run with `--double`.

## Problem explanation

Test problem: `example-problems/50x50x50x1.json`

```json
{
    "dims": 3, "nx": 50, "ny": 50, "nz": 50,
    "substrates_count": 1, "iterations": 1,
    "dx": 20, "dy": 20, "dz": 20,
    "dt": 0.01,
    "diffusion_coefficients": 10000,
    "decay_rates": 0.1,
    "initial_conditions": 1000
}
```

This yields λ = D·dt/dx² = 0.25, giving the Thomas coefficients:
- `b'[i]` (precomputed inverse diagonal): ∈ [0.50, 0.72]
- `c[s]` (off-diagonal): 0.25
- `e[i] = c · b'[i-1]`: ∈ [0.13, 0.18]
- density range during forward sweep: [0, ~1220], annotated as [0, 1300]


## What Was Adapted

### Files changed

| File | Change |
|---|---|
| `src/least_compute_thomas_solver.cpp` | Excluded from `diffuse-taffo` build |
| `taffo-annotations/least_compute_taffo.cpp` | Replaces the above for the TAFFO binary |
| `taffo-annotations/taffo_kernels_v11.cpp` | New file, the annotated kernel, compiled with `taffo` |
| `taffo-annotations/taffo_kernels_v11.h` | New file, declarations |
| `CMakeLists.txt:176–206` | Swaps source file, links kernel object |

`src/least_compute_thomas_solver.h`, all other `src/` files, and the rest of the build
are unchanged.

### What was kept identical

- `precompute_values()` — copied verbatim into `least_compute_taffo.cpp:424–485`
- `initialize()` — copied verbatim (`least_compute_taffo.cpp:487–496`)
- `get_diagonal_layout()` — copied verbatim (`least_compute_taffo.cpp:498–502`)
- `float` path of all `solve_*` methods — uses plain C++ flat-index equivalents of
  the same algorithm; no TAFFO involvement
- OMP structure of `solve()` — same single-region pattern with barriers between dimensions

### What was replaced

The four noarr-based free functions are replaced by flat-index equivalents in
`least_compute_taffo.cpp`. These are generic templates (`generic_x_slice`,
`generic_y_slice_3d`, `generic_z_slice`, etc., lines 42–170) used by the `float` path,
and TAFFO-annotated kernel calls used by the `double` path.

The `double` path for 3D `solve()` calls:
```cpp
taffo_x_slice_v11(dens, bx, cx, ex, ns, nx, M, yz);
taffo_y_row_3d_v11(dens, by, cy, ey, ns, nx, ny, nz, z, x);
taffo_z_row_v11(dens, bz, cz, ez, ns, nx, ny, nz, y, x);
```

## The Annotated Kernel (`taffo_kernels_v11.cpp`)

This file implements the same three-phase Thomas algorithm as the original free functions,
but with two differences:

1. Flat index arithmetic instead of noarr, so the accesses are explicit.
2. A local array with range annotations holds each row's density during the solve.

The annotation strategy applied to each kernel function (example from `taffo_x_slice_v11`,
`taffo_kernels_v11.cpp:61`):

```cpp
// Local array for one row's density — annotated with the value range.
// TAFFO uses this to determine the fixed-point representation.
float local_dens[MAX_NS * MAX_N]
    __attribute__((annotate("scalar(range(0,1300) final)")));
```

The Thomas logic itself is structurally identical to the original. Comparing the forward sweep:

**Original** (`least_compute_thomas_solver.cpp:142–150`, via noarr):
```cpp
dens[yz, i, s] = dens[yz, i, s] + e[i-1, s] * dens[yz, i-1, s];
```

**Annotated kernel** (`taffo_kernels_v11.cpp:69–76`, flat):
```cpp
double d_i   __attribute__((annotate("scalar(range(0,1300) final)")));
double d_prv __attribute__((annotate("scalar(range(0,1300) final)")));
double e_val __attribute__((annotate("scalar(range(0,0.25) final)")));
d_i   = local_dens[s * nx + i];
d_prv = local_dens[s * nx + (i - 1)];
e_val = e[s * nx + (i - 1)];
local_dens[s * nx + i] = (float)(d_i + e_val * d_prv);
```

The range annotations on the scalar temporaries and the local array provide TAFFO with
the value bounds it needs. The copy-out uses `scalar(disabled)` to force conversion back
to double before writing to the caller's `double*` buffer:

```cpp
// taffo_kernels_v11.cpp:98–104
for (int s = 0; s < ns; s++)
    for (int i = 0; i < nx; i++) {
        double out __attribute__((annotate("scalar(disabled)")));
        out = local_dens[s * nx + i];
        dens[s * nx * M + i * M + yz] = out;
    }
```

Each kernel function (x, y-3D, y-2D, z-3D) follows this same pattern. Max sizes are
`MAX_NS = 16` substrates, `MAX_N = 512` elements per dimension.

---

## Performance Results

Benchmark: `--benchmark --double`, 50³ grid, `params.json` (`benchmark_kind: per_dimension`,
5 outer × 10 inner runs, 1 s warmup). Median times.

| Solver | x | y | z | Total | vs GCC |
|---|---|---|---|---|---|
| GCC `lstc --double` (original) | ~115 µs | ~139 µs | ~148 µs | ~402 µs | 1.00× |
| `diffuse-taffo lstc --double` | ~75 µs | ~199 µs | ~65 µs | ~339 µs | 1.19× |


**x and z sweeps** are faster because reduced-precision storage (4 bytes per element vs 8)
halves the working-set size, fitting more data in cache for the sequential/x-strided access
patterns.

**y sweep is slower** because the access stride is `nz` elements = 200 bytes, which exceeds
a cache line (64 bytes). Halving the element size does not reduce cache-miss count at this
stride; the copy-in/copy-out overhead (~30–40 µs) is not recovered.


## Accuracy Results

All RMSE values are deviations from the double-precision baseline: GCC `lstc --double`.
That baseline deviates from the double reference solver by ~1e-13 (fp64 noise only)
and is the ground truth for all comparisons. `diffuse-taffo` is always run with `--double`.

### Single iteration

| Solver | RMSE vs GCC `lstc --double` |
|---|---|
| GCC `lstc --double` (baseline) | ~1e-13 |
| GCC `lstc` (no `--double`, float32) | ~3.4e-05 |
| `diffuse-taffo lstc --double` | ~2.97e-05 |

At a single iteration the annotated solver is slightly more accurate than plain float32,
both deviating from the double baseline by ~3e-05.

### Multi-iteration (full-solve validate, 50³)

RMSE vs GCC `lstc --double` after N full ADI iterations.
Problem: 50³ grid, `init_dens = 1000`, `benchmark_kind: full_solve`.

| Iterations | GCC double (baseline) | GCC float32 | TAFFO |
|---|---|---|---|
| 1   | ~1e-13 | 3.36e-05 | 5.59e-05 |
| 5   | ~1e-13 | 7.33e-05 | 6.36e-05 |
| 10  | ~1e-13 | 2.89e-04 | 3.08e-05 |
| 50  | ~1e-13 | 1.88e-03 | 1.67e-04 |
| 100 | ~1e-13 | 3.56e-03 | 2.01e-04 |
| 500 | ~1e-13 | 1.19e-02 | 1.03e-03 |

Both float32 and TAFFO accumulate error relative to the double baseline, but at different
rates. Float32 grows nearly linearly (N^0.945); TAFFO grows sub-linearly (N^0.851) with a
~10× smaller coefficient. TAFFO becomes more accurate than float32 from around N ≈ 5
iterations onward.

### Grid-size accuracy (50³ / 100³ / 200³, multi-iteration)

All RMSE values vs GCC `lstc --double` on the same grid.

| Grid | Iterations | GCC double (baseline) | GCC float32 | TAFFO |
|---|---|---|---|---|
| 50³  | 10  | ~1e-13 | 2.89e-04 | 3.08e-05 |
| 50³  | 100 | ~1e-13 | 3.56e-03 | 2.01e-04 |
| 100³ | 10  | ~1e-13 | 2.96e-04 | 3.08e-05 |
| 100³ | 100 | ~1e-13 | 3.65e-03 | 2.01e-04 |
| 200³ | 10  | ~1e-13 | 3.01e-04 | 3.08e-05 |
| 200³ | 100 | ~1e-13 | 3.69e-03 | 2.01e-04 |

TAFFO RMSE is grid-size independent, each Thomas row is processed independently, so
grid dimensions affect only element count, not per-row error. Float32 RMSE decreases
slightly with grid size because errors are boundary-concentrated and get diluted by the
larger bulk volume, but both solvers are consistently measured against the double baseline
on the same grid.


## How to Build

### Requirements

- `taffo` on PATH (LLVM-18 frontend, installed separately)
- `cmake` ≥ 3.18
- `clang++-18`
- OpenMP
- LAPACK / BLAS

### Step 1 — Compile the annotated kernel with TAFFO

```bash
mkdir -p scratch

taffo -O3 -std=gnu++20 -fopenmp -march=native \
      -Xdta --maxtotalbits=32 \
      -I taffo-annotations \
      -c taffo-annotations/taffo_kernels_v11.cpp \
      -o scratch/taffo_kernels_v11.o
```
If you want to place the compiled kernels in another folder, you must account to update 'CMakeLists.txt' at line 198.

### Step 2 — Configure and build with CMake

```bash
cmake -B build-clang18 \
      -DCMAKE_CXX_COMPILER=clang++-18 \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build-clang18 --target diffuse-taffo -j$(nproc)
```

### Step 3 — Run

```bash
# Validate correctness vs double reference
./build-clang18/diffuse-taffo \
    --problem example-problems/50x50x50x1.json \
    --alg lstc --params example-problems/params.json \
    --validate --double

# Per-dimension benchmark
./build-clang18/diffuse-taffo \
    --problem example-problems/50x50x50x1.json \
    --alg lstc --params example-problems/params.json \
    --benchmark --double

# Multi-iteration accuracy (full_solve mode)
./build-clang18/diffuse-taffo \
    --problem example-problems/50x50x50x1.json \
    --alg lstc --params example-problems/params_fullsolve.json \
    --validate --double
```

GCC baseline for comparison:
```bash
./build/diffuse \
    --problem example-problems/50x50x50x1.json \
    --alg lstc --params example-problems/params.json \
    --benchmark --double
```

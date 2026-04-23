# TAFFO `lstc` Solver — x86 vs RISC-V Comparison

Problem: `example-problems/50x50x50x1.json` (50³ grid, 1 substrate).  
Benchmark params: `params.json` — `per_dimension`, 5 outer × 10 inner iterations, 1 s warmup.  
Accuracy params: `params_fullsolve.json` — `full_solve` mode, N iterations.  
All RMSE values are vs GCC `lstc --double` on the same platform (double-precision ground truth).

---

## Hardware

| | x86 | RISC-V |
|---|---|---|
| Machine | Intel Core i7-6700K | Banana Pi BPI-F3 |
| ISA | x86-64, AVX2 | RV64GCV (Spacemit X60) |
| Compiler (baseline) | GCC, `-march=native` | clang-18 |
| Compiler (TAFFO) | clang-18 / TAFFO | clang-18 / TAFFO (RISC-V target) |
| TAFFO commit | `1163ab393d` | `1163ab393d` |

---

## Performance

Times in µs (median over 5 outer × 10 inner runs). "Total" = x + y + z sweep.

### RISC-V (Banana Pi BPI-F3, Spacemit X60)

| Solver | x (µs) | y (µs) | z (µs) | Total (µs) | vs GCC double |
|---|---|---|---|---|---|
| baseline `lstc --double` | ~1158 | ~558 | ~1814 | ~3530 | 1.00× |
| baseline `lstc` (float32) | ~1109 | ~552 | ~1742 | ~3403 | 1.04× |
| `diffuse-taffo --double` | ~1368 | ~1436 | ~1146 | ~3950 | 0.89× |

## Accuracy

RMSE vs `lstc --double` (double-precision ground truth) after N full ADI iterations.

| Iterations | double (baseline) | float32 | TAFFO (x86) | TAFFO (RISC-V) |
|---|---|---|---|---|
| 1   | ~1e-13 | 3.36e-05 | 5.59e-05 | 5.59e-05 |
| 5   | ~1e-13 | 7.33e-05 | 6.36e-05 | 6.36e-05 |
| 10  | ~1e-13 | 2.89e-04 | 3.08e-05 | 3.08e-05 |
| 50  | ~1e-13 | 1.88e-03 | 1.67e-04 | 1.67e-04 |
| 100 | ~1e-13 | 3.56e-03 | 2.01e-04 | 2.01e-04 |
| 500 | ~1e-13 | 1.19e-02 | 1.03e-03 | 1.03e-03 |

Results are identical across platforms. TAFFO fixed-point arithmetic is deterministic and hardware-independent: the representation (11 integer + 21 fractional bits) and all operations are fully determined at compile time.

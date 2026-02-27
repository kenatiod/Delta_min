# Delta_min

Number theory research program by Ken Clements (Feb 10, 2026).

## Overview

This program investigates the **pi-complete** property of consecutive integer products `n(n+1)`. For each `n`, it computes:

- **omega** — the number of distinct prime divisors of `n(n+1)`
- **Pidx** — the prime index of the greatest prime divisor of `n(n+1)`
- **delta = Pidx − omega** — the count of "missing" prime divisors between 2 and the greatest prime divisor

When `delta = 0`, the product `n(n+1)` has no missing prime divisors in the sequence from 2 up to its greatest prime factor — i.e., the number is **pi-complete**.

The search is conducted over doubling intervals of powers of 2, tracking the minimum delta, maximum omega, and minimum Pidx found in each interval.

## Algorithm

Because `n` and `n+1` are always coprime:
- `omega(n(n+1)) = omega(n) + omega(n+1)`
- `Pidx(n(n+1)) = max(Pidx(n), Pidx(n+1))`

The loop recycles `omega(n+1)` and `Pidx(n+1)` as the next `omega(n)` and `Pidx(n)` to avoid redundant calculations. Dynamic pruning cuts off factorization early when the result cannot improve the current best delta.

## Parallelism

Uses OpenMP for multi-threaded parallel search with dynamic scheduling. Results from each thread are merged at the end of each interval.

## Build

```bash
gcc -O3 -march=native -fopenmp -o Delta_min Delta_min.c -lm
```

## Usage

```bash
./Delta_min [interval_start] [interval_initial_size] [number_of_intervals] [max_Pidx]
```

**Defaults:** start=1, initial size=1, 40 intervals, max Pidx=30

**Example:**
```bash
./Delta_min 1 1 40 30
```

## Output

For each interval, the program prints:
- Interval start value
- Minimum delta found (and the `n`, Pidx, omega at that minimum)
- Maximum omega found (and the `n`)
- Minimum Pidx found (and the `n`)
- Wall-clock time for the interval

A CSV file `delta_min_data.csv` is also written with the full results.

## Background

The primorial `p_k#` is the product of the first `k` primes. With `max_pidx = 30`, the program can correctly cover `n(n+1)` up to `p_30#`:

```
31,610,054,640,417,607,788,145,206,291,543,662,493,274,686,990
```

## License

MIT License — see [LICENSE](LICENSE).

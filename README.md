# Delta_min_N2

Number theory research program by Ken Clements (Started: February 10, 2026, last revised May 6, 2026).

## Overview

This program investigates the **prime-complete** property of consecutive integer products `n(n+1)`. For each `n`, it computes:

- **omega** — the number of distinct prime divisors of `n(n+1)`
- **iota** — the prime index of the greatest prime divisor of `n(n+1)`
- **delta = iota − omega** — the count of "missing" prime divisors between 2 and the greatest prime divisor

When `delta = 0`, the product `n(n+1)` has no missing prime divisors in the sequence from 2 up to its greatest prime factor — i.e., the number is **prime-complete**.

The search is conducted over doubling intervals of powers of 2, tracking the minimum delta, maximum omega, and minimum iota found in each interval.

Although this README file is written for the program Delta_min_N2.c that tracks the product N_2 = n(n+1), it is also intended for the
corresponding programs Delta_min_N3.c and Delta_min_N4.c that track the products n(n+1)(n+2) and n(n+1)(n+2)(n+3).

## Algorithm

Because `n` and `n+1` are always coprime:
- `omega(n(n+1)) = omega(n) + omega(n+1)`
- `iota(n(n+1)) = max(iota(n), iota(n+1))`

The loop recycles `omega(n+1)` and `iota(n+1)` as the next `omega(n)` and `iota(n)` to avoid redundant calculations. Dynamic pruning cuts off factorization early when the result cannot improve the current best delta.

## Parallelism

Uses OpenMP for multi-threaded parallel search with dynamic scheduling. Results from each thread are merged at the end of each interval.

## Build
For Linux:
```bash
gcc-15 -O3 -march=native -fopenmp -o Delta_min_N2 Delta_min_N2.c -lm
```
For Mac:
```bash
brew install libomp
clang -O3 -std=c11 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -o Delta_min_N2 Delta_min_N2.c -L/opt/homebrew/opt/libomp/lib -lomp -lm
```

## Usage

```bash
./Delta_min_N2 [interval_start] [interval_initial_size] [number_of_intervals] [max_Pidx]
```

**Defaults:** start=1, initial size=1, 40 intervals, max Pidx=100

**Example:**
```bash
./Delta_min_N2 1 1 40 100
```
To sample a single trillion wide interval starting at 100 trillion:
```bash
./Delta_min_N2 100000000000000 1000000000000 1 100
```
A Mac Mini will get through a trillion products in about four hours, 
depending on the max_iota setting. An internal flag, "verbose" can be set in
the source code prior to build to provide detailed output.

## Output

For each interval, the program prints:
- Interval start value
- Minimum delta found (and the `n`, iota, omega at that minimum)
- Maximum omega found (and the `n`)
- Minimum iota found (and the `n`)
- Factoring routine statistics (if verbose)
- Wall-clock time for the interval (if verbose)

A CSV file `delta_min_data.csv` is also written with the full results and a checkpoint file is also written so the program can be restarted at the start of the unfinished interval.

## Background

The primorial $P_r$ is the product of the first $r$ primes. For example, with max_iota = 30, the program can correctly cover $n(n+1)$ up to $P_{30}$=

```
31,610,054,640,417,607,788,145,206,291,543,662,493,274,686,990
```
For any number, $k$, omega($k$) cannot exceed $r$, where $P_r$ is the first primorial greater than or equal to $k$. If we apply this principle to $n(n+1)$, as $n$ increases the product will exceed the values of primorials at about $n$ = the square root of the primorial, because for large numbers $n(n+1)$ is almost exactly $n^2$. This tells us that if $n$ is less than sqrt($P_r$), then omega(n(n+1)) is less than $P_r$.
This limit on maximum omega is useful because if, for intervals between the square roots of primorials, the minimum iota(n(n+1)) is greater than $r$, then it is not possible for any $n(n+1)$ in the interval to be prime-complete. When numbers are large enough such that the minimum iota(n(n+1)) over all subsequent intervals, between primorial square roots, is greater than each corresponding $r$, there is no possibility to find another prime-complete product.

The program can be used to find the list of the last $n$ values for each minimum delta. The default run shows the last $minimum=0$ at $n = 633{,}555$, last $minimum=1$ at $n = 80{,}061{,}344$ and last $minimum=2$ at $n = 1{,}109{,}496{,}723{,}125$. However, no $minimum=3$ appears after $minimum=2$. This is because the $minimum=3$ location is between the $minimum=2$ location and the end of the interval, so it is overshadowed by $minimum=2$. To find the true $minimum=3$ value, the interval must be re-calculated starting right after the $minimum=2$ location and through to the end of the interval. When this is done, the true last $minimum=3$ value is found at $n = 1{,}284{,}729{,}638{,}049$. $Minimum=4$ is at 20,628,591,204,480.

Although Delta_min_N2.c was originally written to iterate through multiple intervals that double in size, single intervals may be run starting at any n value with any length by setting the doubling count to one. For example, to catch that $minimum=3$ value mentioned above, run:
                         ./Delta_min_N2 1109496723126 200000000000 1

See the file "Delta_min_Heuristic_Case.md" for a more detailed explanation, and the file "Delta_min_N2_output.txt" for an output example to $n = 10^{15}$. The Python program "verify_delta_min.py" has been added to the repo so the Delta_min_N2_output.txt file can be automatically checked for correct factorizations producting values for omega and iota.

## Primorial Intervals
The examination of iota is particularly indicative when intervals are set to primorial numbers because of the limits placed on omega there mentioned above. Reaching high primorial intervals is difficult for $N_2$ because iota($N_2$) grows slowly and has large "dips" in early ranges. However for $N_3$ and $N_4$, the iota growth rate is much faster so even a simple Python program on a laptop can reach ranges where iota_min is consistently greater than $r$. Those Python programs and their outputs are included in the N3 and N4 subdirectories.



## License

MIT License — see [LICENSE](LICENSE).

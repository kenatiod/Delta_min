// Delta_min_N4.c

/*
Version 1.0  (derived from Delta_min_N2.c Version 3.0)

This program goes through doubling intervals of n and computes the difference (delta)
between pi(GPF) of the product N_4 = n(n+1)(n+2)(n+3), and its Omega (number of
distinct prime divisors). The prime index of the greatest prime factor, pi(GPF), is
called the "Iota" in this code.

Among any 4 consecutive integers, exactly 2 are even (either n & n+2, or n+1 & n+3),
so the prime factor 2 is always shared and must be subtracted once from the sum of
individual Omega values. Additionally, since n+3 ≡ n (mod 3), both n and n+3 are
divisible by 3 if and only if n ≡ 0 (mod 3); in that case the prime factor 3 is also
shared and requires a second subtraction. No prime p ≥ 5 can divide two of the four
integers (since p > 4 means any two of the four differ by less than p), so no further
corrections are needed.

In summary:
    shared = 1 + (n % 3 == 0 ? 1 : 0)
    Omega(N_4) = Omega(n) + Omega(n+1) + Omega(n+2) + Omega(n+3) - shared

The Iota of N_4 = n(n+1)(n+2)(n+3) is the maximum Iota of the four integers. Computing
those values for n, n+1, n+2, and n+3 individually avoids forming the full product
(which could exceed 256 bits) and keeps all arithmetic in 64-bit integers. Progressing
through successive quadruples, values computed for n+1 through n+3 are reused in the
next quadruple, so each integer is factored only once.

For each doubling interval, the processing looks for the maximum Omega and minimum Iota
for that interval, as well as the minimum delta (Iota - Omega). Once the minimum Iota
grows to exceed the maximum Omega in all successive intervals, no product n(n+1)(n+2)(n+3)
can be prime-complete in those intervals.

At the start of the main routine the variable "verbose" is currently turned off. Set
it to 1 to add statistic details to the printout.

On Linux, compile with:
    gcc -O3 -std=c11 -fopenmp -o Delta_min_N4 Delta_min_N4.c -lm

On Mac, try:
    gcc-15 -O3 -std=c11 -fopenmp -o Delta_min_N4 Delta_min_N4.c -lm

or if needed:
    brew install libomp
    clang -O3 -std=c11 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -o Delta_min_N4 Delta_min_N4.c -L/opt/homebrew/opt/libomp/lib -lomp -lm

You can run it as:
    ./Delta_min_N4                         (starts at 1, runs 40 intervals, uses the first 200 prime divisors)
    ./Delta_min_N4  1000 1000 20 30        (starts at 1000, runs 20 doublings of size starting at 1000, uses 30 primes)
    ./Delta_min_N4  10000000000000 1000000000000  1  300   (starts at 10T, runs 1 interval of size 1T, uses 300 primes)

By Ken Clements, based on Delta_N2_min.c (April 10, 2026). N_4 adaptation May 2026.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_max_threads(void) { return 1; }
static int omp_get_thread_num(void) { return 0; }
static double omp_get_wtime(void) { return 0.0; }
#endif

#define CHECKPOINT_FILE "Delta_min_N4_checkpoint.dat"
#define DEFAULT_BLOCK_N 524288ULL
#define SENTINEL_IOTA 9999999

// Using fewer primes goes faster but gives errors in maxOmega for large numbers
static uint64_t primes[304] = {0,
     2, 3, 5, 7, 11, 13, 17, 19, 23, 29,                            //  10
     31, 37, 41, 43, 47, 53, 59, 61, 67, 71,                        //  20
     73, 79, 83, 89, 97, 101, 103, 107, 109, 113,                   //  30
     127, 131, 137, 139, 149, 151, 157, 163, 167, 173,              //  40
     179, 181, 191, 193, 197, 199, 211, 223, 227, 229,              //  50
     233, 239, 241, 251, 257, 263, 269, 271, 277, 281,              //  60
     283, 293, 307, 311, 313, 317, 331, 337, 347, 349,              //  70
     353, 359, 367, 373, 379, 383, 389, 397, 401, 409,              //  80
     419, 421, 431, 433, 439, 443, 449, 457, 461, 463,              //  90
     467, 479, 487, 491, 499, 503, 509, 521, 523, 541,              // 100
     547, 557, 563, 569, 571, 577, 587, 593, 599, 601,              // 110
     607, 613, 617, 619, 631, 641, 643, 647, 653, 659,              // 120
     661, 673, 677, 683, 691, 701, 709, 719, 727, 733,              // 130
     739, 743, 751, 757, 761, 769, 773, 787, 797, 809,              // 140
     811, 821, 823, 827, 829, 839, 853, 857, 859, 863,              // 150
     877, 881, 883, 887, 907, 911, 919, 929, 937, 941,              // 160
     947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,            // 170
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,    // 180
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,    // 190
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,    // 200
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,    // 210
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,    // 220
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,    // 230
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,    // 240
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,    // 250
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,    // 260
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,    // 270
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,    // 280
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,    // 290
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,    // 300
     1993, 1997, 1999};



static int primes_length = 303; // Only a part of list is used for
        // n(n+1)(n+2)(n+3). Larger list items used for longer consecutive integer products.

// I use this structure to return the Iota and Omega values from the factorization routine.
// If the residue makes it all the way down to 1, the values are exact.
// Even if not exact, it is not possible for the true Iota to be less than
// the value returned (important when looking for minimum Iota over a range).
typedef struct {
    int Iota_small;
    int Omega_small;
    uint64_t residue;
    int exact_Iota;
    int exact_Omega;
} io_result_t;

// Here is the structure for holding the result after a search of an interval.
typedef struct {
    int max_Omega;
    uint64_t mo_n;
    int min_Iota;
    uint64_t mp_n;
    int min_delta;
    int md_Iota;
    int md_Omega;
    uint64_t md_n;
    uint64_t near_max_pairs;
    uint64_t residue_prime_checks;
    uint64_t square_checks;
    uint64_t rho_calls;
    uint64_t full_factor_calls;
} thread_result_t;

// Here is a program flow structure so we can restart after checkpoints.
typedef struct {
    int next_loop_index;
    uint64_t next_n_start;
    uint64_t next_interval_size;
    uint64_t orig_interval_start;
    uint64_t orig_initial_size;
    int orig_intervals;
    int effective_intervals;
    int max_Iota;
    int margin;
    uint64_t block_n;
} checkpoint_t;

// Helper routine to calculate the greatest common denominator.
static inline uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// Helper routine to multiply a times b mod m.
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
#if defined(__SIZEOF_INT128__)
    return (uint64_t)((__uint128_t)a * b % m);
#else
    uint64_t result = 0;
    a %= m;
    while (b) {
        if (b & 1) result = (result + a) % m;
        a = (2 * a) % m;
        b >>= 1;
    }
    return result;
#endif
}

// Helper routine to do modulo exponentiation.
static inline uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp) {
        if (exp & 1) result = mulmod64(result, base, mod);
        base = mulmod64(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// Helper routine that is used in the Miller-Rabin primality test.
static inline int miller_rabin_witness(uint64_t n, uint64_t a, uint64_t d, int r) {
    uint64_t x = powmod64(a, d, n);
    if (x == 1 || x == n - 1) return 1;
    for (int i = 0; i < r - 1; i++) {
        x = mulmod64(x, x, n);
        if (x == n - 1) return 1;
    }
    return 0;
}

// Miller-Rabin primality test for 64 bit integers.
static int is_prime_u64(uint64_t n) {
    if (n < 2) return 0;
    static const uint64_t small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (size_t i = 0; i < sizeof(small_primes)/sizeof(small_primes[0]); i++) {
        if (n == small_primes[i]) return 1;
        if (n % small_primes[i] == 0) return 0;
    }
    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }
    static const uint64_t bases[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (size_t i = 0; i < sizeof(bases)/sizeof(bases[0]); i++) {
        uint64_t a = bases[i];
        if (a >= n) continue;
        if (!miller_rabin_witness(n, a, d, r)) return 0;
    }
    return 1;
}

// Special primality test when the first screening has already been done.
static int is_prime_postsieved_u64(uint64_t n) {

    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }
    static const uint64_t bases[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (size_t i = 0; i < sizeof(bases)/sizeof(bases[0]); i++) {
        uint64_t a = bases[i];
        if (a >= n) continue;
        if (!miller_rabin_witness(n, a, d, r)) return 0;
    }
    return 1;
}



static uint64_t splitmix64_next(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t pollard_rho_single(uint64_t n, uint64_t *seed_state, uint64_t *rho_calls) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;
    (*rho_calls)++;
    while (1) {
        uint64_t c = (splitmix64_next(seed_state) % (n - 1)) + 1;
        uint64_t x = (splitmix64_next(seed_state) % (n - 2)) + 2;
        uint64_t y = x;
        uint64_t d = 1;
        while (d == 1) {
            x = (mulmod64(x, x, n) + c) % n;
            y = (mulmod64(y, y, n) + c) % n;
            y = (mulmod64(y, y, n) + c) % n;
            uint64_t diff = (x > y) ? (x - y) : (y - x);
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
}

// Integer square root routine for 64 bit numbers.
static uint64_t isqrt_u64(uint64_t n) {
    uint64_t x = (uint64_t)sqrt((double)n);
    while ((x + 1) > x && (x + 1) <= UINT64_MAX / (x + 1) && (x + 1) * (x + 1) <= n) x++;
    while (x * x > n) x--;
    return x;
}

/*
This is the central working routine for this program. It is a modified factorization algorithm that
is going to try prime divisors in order, and if they divide n, division continues to remove that
divisor from the residue. The routine counts the divisors that are successful to return the Omega
value, and the returned Iota is the index of the greatest prime that divided n.
A limit is passed to cut off trial divisor attempts at a maximum Iota. If the residue of the trial
divisions goes to 1, the process is stopped and exact Omega and Iota values are returned. However,
if the limit is hit while the residue is still greater than 1, the routine checks to see if the
residue is prime, and if so adds one to the Omega value (which would then be exact) even though the
returned Iota is less than the exact value.
*/
static io_result_t find_io_residue(uint64_t n, int limit_idx) {
    io_result_t out;
    out.Iota_small = 0;
    out.Omega_small = 0;
    out.residue = n;
    out.exact_Iota = 0;
    out.exact_Omega = 0;

    for (int i = 1; i <= limit_idx; i++) {
        uint64_t p = primes[i];
        if (out.residue % p == 0) {
            out.Omega_small++;
            out.Iota_small = i;
            do {
                out.residue /= p;
            } while (out.residue % p == 0);
            if (out.residue == 1) {
                out.exact_Iota = 1;
                out.exact_Omega = 1;
                return out;
            }
        }
    }

    if (out.residue == 1) {
        out.exact_Iota = 1;
        out.exact_Omega = 1;
        return out;
    }

    if (out.residue <= (uint64_t)primes[limit_idx]) {
        for (int i = 1; i <= limit_idx; i++) {
            if (primes[i] == out.residue) {
                out.Omega_small++;
                out.Iota_small = i;
                out.residue = 1;
                out.exact_Iota = 1;
                out.exact_Omega = 1;
                return out;
            }
        }
    }

    return out;
}

// Helper routine to try some more divisors on the residue.
static int upper_bound_extra_factors(uint64_t residue, uint64_t pnext) {
    if (residue <= 1) return 0;
    int cnt = 0;
    while (residue >= pnext) {
        residue /= pnext;
        cnt++;
    }
    return cnt;
}


static int count_distinct_prime_factors_recursive(uint64_t n, uint64_t *seed_state, uint64_t *rho_calls, uint64_t *full_factor_calls) {
    if (n <= 1) return 0;
    if (is_prime_u64(n)) return 1;
    (*full_factor_calls)++;
    uint64_t d = pollard_rho_single(n, seed_state, rho_calls);
    uint64_t a = d;
    uint64_t b = n / d;
    if (a == b) {
        return count_distinct_prime_factors_recursive(a, seed_state, rho_calls, full_factor_calls);
    }
    return count_distinct_prime_factors_recursive(a, seed_state, rho_calls, full_factor_calls) +
           count_distinct_prime_factors_recursive(b, seed_state, rho_calls, full_factor_calls);
}

static int exact_residue_distinct_count(uint64_t residue, uint64_t *seed_state,
                                        uint64_t *prime_checks, uint64_t *square_checks,
                                        uint64_t *rho_calls, uint64_t *full_factor_calls) {
    if (residue <= 1) return 0;

    (*prime_checks)++;
    if (is_prime_postsieved_u64(residue)) return 1;

    (*square_checks)++;
    uint64_t r = isqrt_u64(residue);
    if (r * r == residue) {
        (*prime_checks)++;
        if (is_prime_u64(r)) return 1;
    }

    return count_distinct_prime_factors_recursive(residue, seed_state, rho_calls, full_factor_calls);
}

static char *fmt_num_u64(uint64_t n, char *buf, size_t size) {
    char temp[64];
    snprintf(temp, sizeof(temp), "%llu", (unsigned long long)n);
    int len = (int)strlen(temp);
    int out_idx = 0;
    int first_group = ((len - 1) % 3) + 1;
    for (int i = 0; i < len; i++) {
        if (i > 0 && i >= first_group && ((i - first_group) % 3) == 0) {
            if ((size_t)out_idx < size - 1) buf[out_idx++] = ',';
        }
        if ((size_t)out_idx < size - 1) buf[out_idx++] = temp[i];
    }
    buf[out_idx] = '\0';
    return buf;
}

// Helper routine to format the time field on the printout
static void fmt_time(double sec, char *buf, size_t size) {
    if (sec < 120.0) snprintf(buf, size, "%6.2fs", sec);
    else if (sec < 7200.0) snprintf(buf, size, "%6.2fm", sec / 60.0);
    else if (sec < 172800.0) snprintf(buf, size, "%6.2fh", sec / 3600.0);
    else snprintf(buf, size, "%6.2fd", sec / 86400.0);
}

// Routine to record a checkpoint in the checkpoint file.
static void write_checkpoint(const checkpoint_t *cp) {
    FILE *f = fopen(CHECKPOINT_FILE, "w");
    if (!f) {
        fprintf(stderr, "Warning: could not save checkpoint\n");
        return;
    }
    fprintf(f,
        "orig_interval_start=%llu\n"
        "orig_initial_size=%llu\n"
        "orig_intervals=%d\n"
        "effective_intervals=%d\n"
        "max_Iota=%d\n"
        "margin=%d\n"
        "block_n=%llu\n"
        "next_loop_index=%d\n"
        "next_n_start=%llu\n"
        "next_interval_size=%llu\n",
        (unsigned long long)cp->orig_interval_start,
        (unsigned long long)cp->orig_initial_size,
        cp->orig_intervals,
        cp->effective_intervals,
        cp->max_Iota,
        cp->margin,
        (unsigned long long)cp->block_n,
        cp->next_loop_index,
        (unsigned long long)cp->next_n_start,
        (unsigned long long)cp->next_interval_size);
    fclose(f);
}

// Routine to fetch a checkpoint from the checkpoint file.
static int read_checkpoint(uint64_t orig_interval_start, uint64_t orig_initial_size,
                           int orig_intervals, int max_Iota, int margin, uint64_t block_n,
                           checkpoint_t *cp) {
    FILE *f = fopen(CHECKPOINT_FILE, "r");
    if (!f) return 0;
    memset(cp, 0, sizeof(*cp));
    char line[256];
    int fields = 0;
    while (fgets(line, sizeof(line), f)) {
        unsigned long long ull;
        int iv;
        if      (sscanf(line, "orig_interval_start=%llu", &ull) == 1) { cp->orig_interval_start = (uint64_t)ull; fields++; }
        else if (sscanf(line, "orig_initial_size=%llu", &ull) == 1)   { cp->orig_initial_size = (uint64_t)ull; fields++; }
        else if (sscanf(line, "orig_intervals=%d", &iv) == 1)         { cp->orig_intervals = iv; fields++; }
        else if (sscanf(line, "effective_intervals=%d", &iv) == 1)    { cp->effective_intervals = iv; fields++; }
        else if (sscanf(line, "max_Iota=%d", &iv) == 1)               { cp->max_Iota = iv; fields++; }
        else if (sscanf(line, "margin=%d", &iv) == 1)                 { cp->margin = iv; fields++; }
        else if (sscanf(line, "block_n=%llu", &ull) == 1)             { cp->block_n = (uint64_t)ull; fields++; }
        else if (sscanf(line, "next_loop_index=%d", &iv) == 1)        { cp->next_loop_index = iv; fields++; }
        else if (sscanf(line, "next_n_start=%llu", &ull) == 1)        { cp->next_n_start = (uint64_t)ull; fields++; }
        else if (sscanf(line, "next_interval_size=%llu", &ull) == 1)  { cp->next_interval_size = (uint64_t)ull; fields++; }
    }
    fclose(f);
    if (fields != 10) return 0;
    if (cp->orig_interval_start != orig_interval_start) return 0;
    if (cp->orig_initial_size != orig_initial_size) return 0;
    if (cp->orig_intervals != orig_intervals) return 0;
    if (cp->max_Iota != max_Iota) return 0;
    if (cp->margin != margin) return 0;
    if (cp->block_n != block_n) return 0;
    return 1;
}

/*
Here is the main program routine. We set the default parameters and check to see if the user
changed any of them on the command line. We are going to use as many CPU threads to do the
computation as we can get, so we have to find out how many that is going to be and chop up
the interval blocks to pass them out to the threads. When the threads come back with
their individual results, those results have to be combined to produce the output results
for the interval.
*/
int main(int argc, char *argv[]) {
    uint64_t interval_start = 1;
    uint64_t interval_size = 1;
    int intervals = 40;
    int max_Iota = 200;
    int margin = 3;
    int verbose = 0; // Set this to 1 to get stats in the printout
    uint64_t block_n = DEFAULT_BLOCK_N;

    if (argc > 1) interval_start = strtoull(argv[1], NULL, 10);
    if (argc > 2) interval_size = strtoull(argv[2], NULL, 10);
    if (argc > 3) intervals = atoi(argv[3]);
    if (argc > 4) max_Iota = atoi(argv[4]);
    if (argc > 5) margin = atoi(argv[5]);
    if (argc > 6) block_n = strtoull(argv[6], NULL, 10);
    if (argc > 7) verbose = atoi(argv[5]);

    if (intervals < 1) intervals = 1;
    if (interval_start == 1 && interval_size == 1 && intervals > 1) intervals++;
    if (max_Iota > primes_length) max_Iota = primes_length;
    if (max_Iota < 1) max_Iota = 1;
    if (margin < 0) margin = 0;
    if (block_n < 1024) block_n = 1024;

    const uint64_t pnext = (max_Iota < primes_length) ? primes[max_Iota + 1] : (primes[max_Iota] + 2);

    uint64_t orig_interval_start = interval_start;
    uint64_t orig_initial_size = interval_size;
    int orig_intervals = intervals;

    int resuming = 0;
    int first_interval = 0;
    uint64_t resume_n_stop = 0;
    int effective_intervals = intervals;

    checkpoint_t cp;
    if (read_checkpoint(orig_interval_start, orig_initial_size, orig_intervals, max_Iota, margin, block_n, &cp)) {
        resuming = 1;
        first_interval = cp.next_loop_index;
        resume_n_stop = cp.next_n_start;
        interval_size = cp.next_interval_size;
        effective_intervals = cp.effective_intervals;
    }

    int nthreads = omp_get_max_threads();
    char buf1[64];
    printf("\nDelta_min_N4: Calculating the smallest difference between the prime number index of the greatest prime divisor\n");
    if (intervals == 1) {
        printf("     of N_4 = n(n+1)(n+2)(n+3) and the count of prime divisors of that consecutive integer product, over an interval.\n");
    } else {
        printf("     of N_4 = n(n+1)(n+2)(n+3) and the count of prime divisors of that consecutive integer product, over doubling intervals.\n");
    }
    printf("====================================================================================================================\n\n");

    if (resuming) {
        printf("*** Resuming from checkpoint: starting at loop interval %d, n = %llu ***\n\n",
               first_interval, (unsigned long long)resume_n_stop);
    }
    printf("Tracking: max(Omega) if full factorization, exact min(Iota), and exact min(Iota-Omega) per interval.\n");
    printf("Omega(N_4) = Omega(n)+Omega(n+1)+Omega(n+2)+Omega(n+3) - 1 (for shared 2) - (1 if n%%3==0, for shared 3).\n");
    int print_intervals = orig_intervals;
    if (print_intervals > 1) print_intervals--;
    printf("Intervals: %d, initial size %s\n", print_intervals,
           fmt_num_u64(orig_initial_size, buf1, sizeof(buf1)));
    printf("Threads: %d, max prime index: %d (p_%d = %llu), margin: %d, block_n: %llu\n\n",
           nthreads, max_Iota, max_Iota, (unsigned long long)primes[max_Iota],
           margin, (unsigned long long)block_n);
    if (verbose) {
        printf("            Interval Start  minDelta = Iota - Omega      n@minDelta  maxOmega      n@maxOmega  minIota        n@minIota    nearPairs  primeChk   rhoCalls  deepFact   time\n");
        printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    } else {
        printf("            Interval Start  minDelta = Iota - Omega      n@minDelta  maxOmega      n@maxOmega  minIota        n@minIota \n");
        printf("---------------------------------------------------------------------------------------------------------------------------\n");
    }
    fflush(stdout);

    uint64_t n_stop = 0;
    if (resuming) {
        intervals = effective_intervals;
        n_stop = resume_n_stop;
    } else {
        // Special printout for the case of an interval of 1 number starting at 1.
        // N_4(1) = 1*2*3*4 = 24 = 2^3 * 3.  Distinct primes: {2,3}.
        // omega(1)+omega(2)+omega(3)+omega(4) = 0+1+1+1 = 3.
        // n=1: not divisible by 3, so shared=1 (only factor 2).
        // Omega = 3 - 1 = 2. Iota = max(0,1,2,1) = pi(3) = 2. delta = 0.
        if (orig_interval_start == 1 && orig_initial_size == 1) {
            if (verbose) {
                printf("                         1         0     2     2                  1      2                  1      2                  1           0         0         0         0    0.00s\n");
            } else {
                printf("                         1         0     2     2                  1      2                  1      2                  1  \n");
            }
            first_interval = 1;
            interval_start = 2;
            interval_size = 2;
            intervals--;
            n_stop = 2;
        }
        effective_intervals = intervals;
    }

    FILE *csv = NULL;
    if (resuming) {
        csv = fopen("Delta_min_N4_data.csv", "a");
    } else {
        csv = fopen("Delta_min_N4_data.csv", "w");
        if (csv) {
            fprintf(csv, "interval,n_start,n_stop,min_delta,md_Iota,md_Omega,md_n,max_Omega,mo_n,min_Iota,mp_n,margin,near_max_pairs,residue_prime_checks,square_checks,rho_calls,full_factor_calls\n");
            if (orig_interval_start == 1 && orig_initial_size == 1) {
                // N_4(1)=24: Omega=2, Iota=2, delta=0
                fprintf(csv, "0,1,1,0,2,2,1,2,1,2,1,%d,0,0,0,0,0\n", margin);
            }
        }
    }

    for (int interval = first_interval; interval < intervals; interval++) {
        uint64_t n_start = (interval == 0) ? interval_start : n_stop;
        n_stop = n_start + interval_size;
        interval_size *= 2;

        double t0 = omp_get_wtime();
        thread_result_t *thr = (thread_result_t *)calloc((size_t)nthreads, sizeof(thread_result_t));
        if (!thr) {
            fprintf(stderr, "Allocation failure for thread results\n");
            return 1;
        }
        for (int t = 0; t < nthreads; t++) {
            thr[t].max_Omega = -1;
            thr[t].min_Iota = SENTINEL_IOTA;
            thr[t].min_delta = SENTINEL_IOTA;
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int tid = omp_get_thread_num();
            thread_result_t *my = &thr[tid];
            uint64_t seed_state = 0x123456789abcdef0ULL ^ ((uint64_t)(tid + 1) * 0x9E3779B97F4A7C15ULL);

            // We need n, n+1, n+2, n+3 for each quadruple, so we allocate three extra slots.
            uint64_t total_numbers = (n_stop + 3) - n_start + 1;
            uint64_t num_blocks = (total_numbers + block_n - 1) / block_n;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (uint64_t blk = 0; blk < num_blocks; blk++) {
                uint64_t block_start = n_start + blk * block_n;
                // Need three extra slots beyond block_n to cover the n+3 element of the last quadruple.
                uint64_t block_end_excl = block_start + block_n + 3;
                if (block_end_excl > n_stop + 4) block_end_excl = n_stop + 4;
                size_t len = (size_t)(block_end_excl - block_start);
                io_result_t *arr = (io_result_t *)malloc(len * sizeof(io_result_t));
                if (!arr) {
                    fprintf(stderr, "Allocation failure for block array\n");
                    exit(1);
                }

                for (size_t i = 0; i < len; i++) {
                    arr[i] = find_io_residue(block_start + (uint64_t)i, max_Iota);
                }

                // Each quadruple is (arr[i], arr[i+1], arr[i+2], arr[i+3]) for n, n+1, n+2, n+3.
                size_t quad_count = (len >= 3) ? len - 3 : 0;
                for (size_t i = 0; i < quad_count; i++) {
                    uint64_t n = block_start + (uint64_t)i;
                    io_result_t *a = &arr[i];       // n
                    io_result_t *b = &arr[i + 1];   // n+1
                    io_result_t *c = &arr[i + 2];   // n+2
                    io_result_t *d = &arr[i + 3];   // n+3

                    // Shared prime correction:
                    //   Factor 2 is always shared (exactly 2 of the 4 are even).
                    //   Factor 3 is shared when n ≡ 0 (mod 3), because then n and n+3 are
                    //   both divisible by 3.  No prime ≥ 5 can be shared among the four.
                    int shared = 1 + (n % 3 == 0 ? 1 : 0);
                    int Omega_lb = a->Omega_small + b->Omega_small + c->Omega_small + d->Omega_small - shared;

                    // Iota: the maximum prime index among the four integers.
                    if (a->exact_Iota && b->exact_Iota && c->exact_Iota && d->exact_Iota) {
                        int Iota = a->Iota_small;
                        if (b->Iota_small > Iota) Iota = b->Iota_small;
                        if (c->Iota_small > Iota) Iota = c->Iota_small;
                        if (d->Iota_small > Iota) Iota = d->Iota_small;

                        if (Iota < my->min_Iota || (Iota == my->min_Iota && n > my->mp_n)) {
                            my->min_Iota = Iota;
                            my->mp_n = n;
                        }
                        if (a->exact_Omega && b->exact_Omega && c->exact_Omega && d->exact_Omega) {
                            int delta = Iota - Omega_lb;
                            if (delta < my->min_delta || (delta == my->min_delta && n > my->md_n)) {
                                my->min_delta = delta;
                                my->md_Iota = Iota;
                                my->md_Omega = Omega_lb;
                                my->md_n = n;
                            }
                        }
                    }

                    if (a->exact_Omega && b->exact_Omega && c->exact_Omega && d->exact_Omega) {
                        if (Omega_lb > my->max_Omega || (Omega_lb == my->max_Omega && n > my->mo_n)) {
                            my->max_Omega = Omega_lb;
                            my->mo_n = n;
                        }
                        continue;
                    }

                    if (Omega_lb < my->max_Omega - margin) {
                        continue;
                    }

                    my->near_max_pairs++;
                    int ub = Omega_lb
                             + upper_bound_extra_factors(a->residue, pnext)
                             + upper_bound_extra_factors(b->residue, pnext)
                             + upper_bound_extra_factors(c->residue, pnext)
                             + upper_bound_extra_factors(d->residue, pnext);
                    if (ub < my->max_Omega) {
                        continue;
                    }

                    int Omega_exact = Omega_lb;
                    if (!a->exact_Omega) {
                        Omega_exact += exact_residue_distinct_count(a->residue, &seed_state,
                                                                    &my->residue_prime_checks,
                                                                    &my->square_checks,
                                                                    &my->rho_calls,
                                                                    &my->full_factor_calls);
                    }
                    if (!b->exact_Omega) {
                        Omega_exact += exact_residue_distinct_count(b->residue, &seed_state,
                                                                    &my->residue_prime_checks,
                                                                    &my->square_checks,
                                                                    &my->rho_calls,
                                                                    &my->full_factor_calls);
                    }
                    if (!c->exact_Omega) {
                        Omega_exact += exact_residue_distinct_count(c->residue, &seed_state,
                                                                    &my->residue_prime_checks,
                                                                    &my->square_checks,
                                                                    &my->rho_calls,
                                                                    &my->full_factor_calls);
                    }
                    if (!d->exact_Omega) {
                        Omega_exact += exact_residue_distinct_count(d->residue, &seed_state,
                                                                    &my->residue_prime_checks,
                                                                    &my->square_checks,
                                                                    &my->rho_calls,
                                                                    &my->full_factor_calls);
                    }

                    if (Omega_exact > my->max_Omega || (Omega_exact == my->max_Omega && n > my->mo_n)) {
                        my->max_Omega = Omega_exact;
                        my->mo_n = n;
                    }
                }

                free(arr);
            }
        }

        thread_result_t best;
        memset(&best, 0, sizeof(best));
        best.max_Omega = -1;
        best.min_Iota = SENTINEL_IOTA;
        best.min_delta = SENTINEL_IOTA;
        for (int t = 0; t < nthreads; t++) {
            if (thr[t].max_Omega > best.max_Omega || (thr[t].max_Omega == best.max_Omega && thr[t].mo_n > best.mo_n)) {
                best.max_Omega = thr[t].max_Omega;
                best.mo_n = thr[t].mo_n;
            }
            if (thr[t].min_Iota < best.min_Iota || (thr[t].min_Iota == best.min_Iota && thr[t].mp_n > best.mp_n)) {
                best.min_Iota = thr[t].min_Iota;
                best.mp_n = thr[t].mp_n;
            }
            if (thr[t].min_delta < best.min_delta || (thr[t].min_delta == best.min_delta && thr[t].md_n > best.md_n)) {
                best.min_delta = thr[t].min_delta;
                best.md_Iota = thr[t].md_Iota;
                best.md_Omega = thr[t].md_Omega;
                best.md_n = thr[t].md_n;
            }
            best.near_max_pairs += thr[t].near_max_pairs;
            best.residue_prime_checks += thr[t].residue_prime_checks;
            best.square_checks += thr[t].square_checks;
            best.rho_calls += thr[t].rho_calls;
            best.full_factor_calls += thr[t].full_factor_calls;
        }
        free(thr);

        double dt = omp_get_wtime() - t0;
        char start_str[64], md_str[64], mo_str[64], mp_str[64], time_str[32];
        fmt_num_u64(n_start, start_str, sizeof(start_str));
        fmt_num_u64(best.md_n, md_str, sizeof(md_str));
        fmt_num_u64(best.mo_n, mo_str, sizeof(mo_str));
        fmt_num_u64(best.mp_n, mp_str, sizeof(mp_str));
        fmt_time(dt, time_str, sizeof(time_str));
        if (verbose) {
            printf("%26s %9d %5d %5d %18s %6d %18s %6d %18s %11llu %9llu %9llu %9llu %8s\n",
               start_str,
               best.min_delta == SENTINEL_IOTA ? 9999 : best.min_delta,
               best.md_Iota,
               best.md_Omega,
               md_str,
               best.max_Omega,
               mo_str,
               best.min_Iota == SENTINEL_IOTA ? 9999 : best.min_Iota,
               mp_str,
               (unsigned long long)best.near_max_pairs,
               (unsigned long long)best.residue_prime_checks,
               (unsigned long long)best.rho_calls,
               (unsigned long long)best.full_factor_calls,
               time_str);
         } else {
            printf("%26s %9d %5d %5d %18s %6d %18s %6d %18s\n",
               start_str,
               best.min_delta == SENTINEL_IOTA ? 9999 : best.min_delta,
               best.md_Iota,
               best.md_Omega,
               md_str,
               best.max_Omega,
               mo_str,
               best.min_Iota == SENTINEL_IOTA ? 9999 : best.min_Iota,
               mp_str);
         }
        fflush(stdout);

        if (csv) {
            fprintf(csv, "%d,%llu,%llu,%d,%d,%d,%llu,%d,%llu,%d,%llu,%d,%llu,%llu,%llu,%llu,%llu\n",
                    interval,
                    (unsigned long long)n_start,
                    (unsigned long long)(n_stop - 1),
                    best.min_delta == SENTINEL_IOTA ? 9999 : best.min_delta,
                    best.md_Iota,
                    best.md_Omega,
                    (unsigned long long)best.md_n,
                    best.max_Omega,
                    (unsigned long long)best.mo_n,
                    best.min_Iota,
                    (unsigned long long)best.mp_n,
                    margin,
                    (unsigned long long)best.near_max_pairs,
                    (unsigned long long)best.residue_prime_checks,
                    (unsigned long long)best.square_checks,
                    (unsigned long long)best.rho_calls,
                    (unsigned long long)best.full_factor_calls);
            fflush(csv);
        }

        checkpoint_t outcp;
        outcp.next_loop_index = interval + 1;
        outcp.next_n_start = n_stop;
        outcp.next_interval_size = interval_size;
        outcp.orig_interval_start = orig_interval_start;
        outcp.orig_initial_size = orig_initial_size;
        outcp.orig_intervals = orig_intervals;
        outcp.effective_intervals = effective_intervals;
        outcp.max_Iota = max_Iota;
        outcp.margin = margin;
        outcp.block_n = block_n;
        write_checkpoint(&outcp);
    }

    if (csv) fclose(csv);
    remove(CHECKPOINT_FILE);
    printf("\nDone.\n");
    return 0;
}

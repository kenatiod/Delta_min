/*
 * Delta_min.c
 *
 * This program checks the product n(n+1) over intervals of powers of 2. The number
 * of prime divisors of the product (omega) is calculated for each n, together with
 * the prime number index of the greatest prime divisor (Pidx) of the product. The
 * significance of these two values is that delta = Pidx - omega gives the count of
 * missing prime divisors so when those are equal, there are no missiong prime divisors
 * in the set from 2 all the way to the greatest prime divisor, i.e. the number is
 * pi-complete.
 *
 * The program uses multiple parallel CPU threads, and each thread computes the
 * omega(n(n+1)) by adding omega(n)+omega(n+1) because n and (n+1) are coprime. To
 * compute Pidx(n(n+1)) it is simply the greater value of Pidx(n) v. Pidx(n+1).
 * In order to speed the processing, the loop recycles the values of Omega(n+1)
 * and Pidx(n+1) so as not to repeate the calculation.
 *
 * The actual calculation of omega and Pidx are limited to values that could
 * impact the recording of the minimum delta. This is valid because each interval
 * is limited to n(n+1) less than or equal to the primorial number p_k#, where
 * k is the maximum omega. For example, if the routine can properly calculate
 * omega to 30, it can cover n(n+1) up to p_30# which is:
 * 31,610,054,640,417,607,788,145,206,291,543,662,493,274,686,990
 *
 * Alternate runs can be produced by command line parameters:
 *   ./Delta.min interval_start interval_initial_size number_of_intervals max_Pidx
 *
 * Checkpoint/resume: the program saves delta_min_checkpoint.dat after each
 * completed interval. If interrupted, restart with the same arguments and it
 * will resume automatically from where it left off. The checkpoint file is
 * removed on successful completion.
 *
 * Compile: gcc-15 -O3 -march=native -fopenmp -o Delta_min Delta_min.c -lm
 *
 * By Ken Clements Feb 10, 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <math.h>

#define MAX_PRIMES 256
#define CHECKPOINT_FILE "delta_min_checkpoint.dat"

static uint64_t primes[MAX_PRIMES + 1];
static int num_primes = 0;

/* Shared pruning variable for the current interval */
volatile int g_min_delta = 999999;

static int is_prime_test(uint64_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (uint64_t d = 5; d * d <= n; d += 6) {
        if (n % d == 0 || n % (d + 2) == 0) return 0;
    }
    return 1;
}

static void build_prime_table(int count) {
    num_primes = 0;
    primes[0] = 0;
    uint64_t c = 2;
    while (num_primes < count) {
        if (is_prime_test(c)) {
            num_primes++;
            primes[num_primes] = c;
        }
        c++;
    }
}

/*
 * Compute Pidx and omega for a single integer.
 * Returns 1 if successful (n is smooth within limit), 0 otherwise.
 */
static inline int find_pidx_omega(uint64_t n, int limit_idx,
                                  int *pidx_out, int *omega_out) {
    int pidx = 0, omega = 0;

    for (int i = 1; i <= limit_idx; i++) {
        uint64_t p = primes[i];
        if (p > n) break;

        if (n % p == 0) {
            omega++;
            pidx = i;
            do { n /= p; } while (n % p == 0);

            if (n == 1) {
                *pidx_out = pidx;
                *omega_out = omega;
                return 1;
            }
        }
    }
    return 0;
}

/* Format number with commas into a provided buffer */
static char* fmt_num(uint64_t n, char *buf, size_t size) {
    char temp[64];
    snprintf(temp, sizeof(temp), "%lu", (unsigned long)n);
    int len = (int)strlen(temp);
    int out_idx = 0;
    int first_group = ((len - 1) % 3) + 1;

    for (int i = 0; i < len; i++) {
        if (i > 0 && ((i - first_group) % 3) == 0 && (i >= first_group)) {
            if (out_idx < size - 1) buf[out_idx++] = ',';
        }
        if (out_idx < size - 1) buf[out_idx++] = temp[i];
    }
    buf[out_idx] = '\0';
    return buf;
}

/* Format time into best units */
static void fmt_time(double sec, char *buf, size_t size) {
    if (sec < 120.0) {
        snprintf(buf, size, "%6.2fs", sec);
    } else if (sec < 7200.0) { // < 2 hours
        snprintf(buf, size, "%6.2fm", sec / 60.0);
    } else if (sec < 172800.0) { // < 48 hours (2 days)
        snprintf(buf, size, "%6.2fh", sec / 3600.0);
    } else {
        snprintf(buf, size, "%6.2fd", sec / 86400.0);
    }
}

/* ---- Checkpoint support ---- */

/*
 * Save run state to CHECKPOINT_FILE so the program can resume after interruption.
 * Called after every completed interval.
 */
static void write_checkpoint(int next_loop_index, uint64_t next_n_start,
                              uint64_t next_interval_size,
                              uint64_t orig_interval_start, uint64_t orig_initial_size,
                              int orig_intervals, int effective_intervals, int max_pidx) {
    FILE *f = fopen(CHECKPOINT_FILE, "w");
    if (!f) { fprintf(stderr, "Warning: could not save checkpoint\n"); return; }
    fprintf(f,
        "orig_interval_start=%llu\n"
        "orig_initial_size=%llu\n"
        "orig_intervals=%d\n"
        "effective_intervals=%d\n"
        "max_pidx=%d\n"
        "next_loop_index=%d\n"
        "next_n_start=%llu\n"
        "next_interval_size=%llu\n",
        (unsigned long long)orig_interval_start,
        (unsigned long long)orig_initial_size,
        orig_intervals,
        effective_intervals,
        max_pidx,
        next_loop_index,
        (unsigned long long)next_n_start,
        (unsigned long long)next_interval_size);
    fclose(f);
}

/*
 * Try to load a checkpoint matching the current run parameters.
 * Returns 1 and populates output args on success, 0 otherwise.
 */
static int read_checkpoint(uint64_t orig_interval_start, uint64_t orig_initial_size,
                            int orig_intervals, int max_pidx,
                            int *out_loop_index, uint64_t *out_n_start,
                            uint64_t *out_interval_size, int *out_effective_intervals) {
    FILE *f = fopen(CHECKPOINT_FILE, "r");
    if (!f) return 0;

    unsigned long long rd_orig_start = 0, rd_orig_size = 0;
    unsigned long long rd_next_n = 0, rd_next_size = 0;
    int rd_orig_ivals = 0, rd_eff_ivals = 0, rd_max_pidx = 0, rd_next_idx = 0;
    int fields = 0;
    char line[256];

    while (fgets(line, sizeof(line), f)) {
        unsigned long long ull; int iv;
        if      (sscanf(line, "orig_interval_start=%llu", &ull) == 1) { rd_orig_start  = ull; fields++; }
        else if (sscanf(line, "orig_initial_size=%llu",   &ull) == 1) { rd_orig_size   = ull; fields++; }
        else if (sscanf(line, "orig_intervals=%d",        &iv)  == 1) { rd_orig_ivals  = iv;  fields++; }
        else if (sscanf(line, "effective_intervals=%d",   &iv)  == 1) { rd_eff_ivals   = iv;  fields++; }
        else if (sscanf(line, "max_pidx=%d",              &iv)  == 1) { rd_max_pidx    = iv;  fields++; }
        else if (sscanf(line, "next_loop_index=%d",       &iv)  == 1) { rd_next_idx    = iv;  fields++; }
        else if (sscanf(line, "next_n_start=%llu",        &ull) == 1) { rd_next_n      = ull; fields++; }
        else if (sscanf(line, "next_interval_size=%llu",  &ull) == 1) { rd_next_size   = ull; fields++; }
    }
    fclose(f);

    if (fields != 8) return 0;

    /* Validate that the checkpoint matches the current run parameters */
    if (rd_orig_start != (unsigned long long)orig_interval_start) return 0;
    if (rd_orig_size  != (unsigned long long)orig_initial_size)   return 0;
    if (rd_orig_ivals != orig_intervals)                          return 0;
    if (rd_max_pidx   != max_pidx)                                return 0;

    *out_loop_index          = rd_next_idx;
    *out_n_start             = (uint64_t)rd_next_n;
    *out_interval_size       = (uint64_t)rd_next_size;
    *out_effective_intervals = rd_eff_ivals;
    return 1;
}

/* ---- Per-interval max_pidx cap via primorial lookup ---- */

/*
 * The k-th primorial p_k# is the product of the first k primes.
 * For a single integer n, omega(n) <= k only when p_k# <= n.
 * Therefore, for any n in [n_start, n_stop], the maximum achievable
 * omega is the largest k where p_k# <= n_stop.  Checking primes beyond
 * that index wastes time on factorizations that cannot contribute to
 * max_omega, min_pidx, or min_delta for this interval.
 *
 * p_16# overflows uint64_t, so the table covers k = 1..15 and the
 * bound is naturally capped at 15 for any uint64_t value of n_stop.
 */
static const uint64_t primorials[16] = {
    0ULL,                     /* index 0 unused */
    2ULL,                     /* p_1#  = 2 */
    6ULL,                     /* p_2#  = 2·3 */
    30ULL,                    /* p_3#  = 2·3·5 */
    210ULL,                   /* p_4#  = 2·3·5·7 */
    2310ULL,                  /* p_5#  = ····11 */
    30030ULL,                 /* p_6#  = ····13 */
    510510ULL,                /* p_7#  = ····17 */
    9699690ULL,               /* p_8#  = ····19 */
    223092870ULL,             /* p_9#  = ····23 */
    6469693230ULL,            /* p_10# = ····29 */
    200560490130ULL,          /* p_11# = ····31 */
    7420738134810ULL,         /* p_12# = ····37 */
    304250263527210ULL,       /* p_13# = ····41 */
    13082761331670030ULL,     /* p_14# = ····43 */
    614889782588491410ULL,    /* p_15# = ····47 */
};

/* Return the largest k (1..15) such that primorials[k] <= n, or 1 if n < 2. */
static int omega_bound(uint64_t n) {
    int k = 1;
    for (int i = 2; i <= 15; i++) {
        if (primorials[i] <= n) k = i;
        else break;
    }
    return k;
}

/* Per-thread result for reduction */
typedef struct {
    int min_delta;
    int md_pidx;
    int md_omega;
    uint64_t md_n;

    int max_omega;
    uint64_t mo_n;

    int min_pidx;
    uint64_t mp_n;
} interval_result_t;

int main(int argc, char *argv[]) {
    uint64_t interval_start = 1;
    uint64_t interval_size = 1;
    int intervals = 40;
    int max_pidx = 30;

    if (argc > 1) interval_start = strtoull(argv[1], NULL, 10);
    if (argc > 2) interval_size = strtoull(argv[2], NULL, 10);
    if (argc > 3) intervals = atoi(argv[3]);
    if (argc > 4) max_pidx = atoi(argv[4]);

    /* Save original command-line params for checkpoint validation */
    uint64_t orig_interval_start = interval_start;
    uint64_t orig_initial_size   = interval_size;
    int      orig_intervals      = intervals;

    int table_size = max_pidx + 50;
    if (table_size > MAX_PRIMES) table_size = MAX_PRIMES;
    build_prime_table(table_size);

    if (max_pidx > num_primes) max_pidx = num_primes;

    /* Check for a checkpoint from a previous interrupted run */
    int resuming = 0;
    int start_interval = 0;
    uint64_t resume_n_stop = 0;
    int effective_intervals = intervals;

    {
        int ci = 0, ce = intervals;
        uint64_t cn = 0, cs = interval_size;
        if (read_checkpoint(orig_interval_start, orig_initial_size, orig_intervals, max_pidx,
                            &ci, &cn, &cs, &ce)) {
            resuming            = 1;
            start_interval      = ci;
            resume_n_stop       = cn;
            interval_size       = cs;
            effective_intervals = ce;
        }
    }

    int nthreads = omp_get_max_threads();
    char buf1[64];

    printf("\nDelta_min: Formatted Output & Optimized Logic\n");
    printf("=================================================\n\n");
    if (resuming) {
        printf("*** Resuming from checkpoint: starting at loop interval %d, n = %llu ***\n\n",
               start_interval, (unsigned long long)resume_n_stop);
    }
    printf("Tracking: min(Pidx-omega), max(omega), min(Pidx) per interval.\n");
    printf("Intervals: %d doubling, initial size %s\n", orig_intervals,
           fmt_num(orig_initial_size, buf1, sizeof(buf1)));
    printf("Threads: %d, max prime index: %d (p_%d = %lu)\n\n",
           nthreads, max_pidx, max_pidx, (unsigned long)primes[max_pidx]);

    printf("            Interval Start  minDelta = Pidx - omega      n@minDelta  maxOmmega     n@maxOmega  minPidx        n@minPidx     time\n");
    printf("--------------------------------------------------------------------------------------------------------------------------------------\n");

    if (resuming) {
        /* Restore the effective interval count (already adjusted for n=1 if applicable) */
        intervals = effective_intervals;
    } else {
        /* Fresh run: handle the n=1 special case */
        if (interval_start == 1 && interval_size == 1) {
            printf("                         1         0     0     0                  0      0                  1      0                  1    0.00s\n");
            interval_start = 2;
            interval_size = 2;
            intervals--;
        }
        effective_intervals = intervals;
    }
    fflush(stdout);

    uint64_t n_start, n_stop;
    /* Pre-set n_stop so the loop's "n_start = n_stop" branch works correctly on resume */
    n_stop = resuming ? resume_n_stop : 0;

    /* Open CSV: append on resume (skip header), write fresh otherwise */
    FILE *csv;
    if (resuming) {
        csv = fopen("delta_min_data.csv", "a");
    } else {
        csv = fopen("delta_min_data.csv", "w");
        if (csv) {
            fprintf(csv, "interval,n_start,n_stop,min_delta,md_pidx,md_omega,md_n,max_omega,mo_n,min_pidx,mp_n\n");
        }
    }

    for (int interval = start_interval; interval < intervals; interval++) {
        if (interval == 0) {
            n_start = interval_start;
        } else {
            n_start = n_stop;
        }
        n_stop = n_start + interval_size;
        interval_size *= 2;

        /* Tighten the prime search limit to what is actually achievable for
         * individual n <= n_stop.  omega(n) cannot exceed k where p_k# <= n_stop,
         * so there is no point checking primes beyond index k.  This replaces
         * max_pidx as the hard ceiling for this interval; g_min_delta+22 pruning
         * still tightens the limit further on a per-n basis as the interval runs. */
        int eff_max_pidx = omega_bound(n_stop);
        if (eff_max_pidx > max_pidx) eff_max_pidx = max_pidx;

        double t0 = omp_get_wtime();
        g_min_delta = 999999;

        interval_result_t *thr_res = calloc(nthreads, sizeof(interval_result_t));
        if (!thr_res) exit(1);

        for (int t = 0; t < nthreads; t++) {
            thr_res[t].min_delta = 999999999;
            thr_res[t].max_omega = -1;
            thr_res[t].min_pidx = 999999999;
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            interval_result_t *my = &thr_res[tid];


            int p1 = 0, o1 = 0;
            if (!find_pidx_omega(n_start, eff_max_pidx, &p1, &o1)){
                p1 = 9999999;
                o1 = 0;
            }
            #pragma omp for schedule(dynamic, 65536)
            for (uint64_t n = n_start+1; n <= n_stop; n++) {

                // Dynamic Pruning
                int current_best = g_min_delta;
                int limit = current_best + 22;
                if (limit > eff_max_pidx) limit = eff_max_pidx;

                int p2 = 0, o2 = 0;
                if (find_pidx_omega(n, limit, &p2, &o2)) {

                    int pidx = (p1 > p2) ? p1 : p2;
                    int omega = o1 + o2;
                    int delta = pidx - omega;

                    if (delta <= my->min_delta) {
                        my->min_delta = delta;
                        my->md_pidx = pidx;
                        my->md_omega = omega;
                        my->md_n = n-1;

                        if (delta < g_min_delta) {
                            #pragma omp atomic write
                            g_min_delta = delta;
                        }
                    }

                    if (omega >= my->max_omega) {
                        my->max_omega = omega;
                        my->mo_n = n-1;
                    }

                    if (pidx <= my->min_pidx) {
                        my->min_pidx = pidx;
                        my->mp_n = n-1;
                    }
                    p1 = p2;
                    o1 = o2;
                } else {
                    p1 = 9999999;
                    o1 = 1;
                }
            }
        }

        /* Merge results */
        interval_result_t best;
        best.min_delta = 999999999;
        best.max_omega = -1;
        best.min_pidx = 999999999;

        for (int t = 0; t < nthreads; t++) {
            if (thr_res[t].min_delta < best.min_delta ||
               (thr_res[t].min_delta == best.min_delta && thr_res[t].md_n > best.md_n)) {
                best.min_delta = thr_res[t].min_delta;
                best.md_pidx = thr_res[t].md_pidx;
                best.md_omega = thr_res[t].md_omega;
                best.md_n = thr_res[t].md_n;
            }
            if (thr_res[t].max_omega > best.max_omega ||
               (thr_res[t].max_omega == best.max_omega && thr_res[t].mo_n > best.mo_n)) {
                best.max_omega = thr_res[t].max_omega;
                best.mo_n = thr_res[t].mo_n;
            }
            if (thr_res[t].min_pidx < best.min_pidx ||
               (thr_res[t].min_pidx == best.min_pidx && thr_res[t].mp_n > best.mp_n)) {
                best.min_pidx = thr_res[t].min_pidx;
                best.mp_n = thr_res[t].mp_n;
            }
        }
        free(thr_res);

        double dt = omp_get_wtime() - t0;

        char start_str[64], md_n_str[64], mo_n_str[64], mp_n_str[64], time_str[32];
        fmt_num(n_start, start_str, sizeof(start_str));
        fmt_num(best.md_n, md_n_str, sizeof(md_n_str));
        fmt_num(best.mo_n, mo_n_str, sizeof(mo_n_str));
        fmt_num(best.mp_n, mp_n_str, sizeof(mp_n_str));
        fmt_time(dt, time_str, sizeof(time_str));

        printf("%26s %9d %5d %5d %18s %6d %18s %6d %18s %8s\n",
               start_str,
               best.min_delta == 999999999 ? 9999 : best.min_delta,
               best.md_pidx, best.md_omega, md_n_str,
               best.max_omega, mo_n_str,
               best.min_pidx == 999999999 ? 9999 : best.min_pidx,
               mp_n_str,
               time_str);
        fflush(stdout);

        if (csv) {
            fprintf(csv, "%d,%lu,%lu,%d,%d,%d,%lu,%d,%lu,%d,%lu\n",
                    interval, (unsigned long)n_start, (unsigned long)n_stop,
                    best.min_delta, best.md_pidx, best.md_omega, (unsigned long)best.md_n,
                    best.max_omega, (unsigned long)best.mo_n,
                    best.min_pidx, (unsigned long)best.mp_n);
            fflush(csv);
        }

        /* Save checkpoint so we can resume if interrupted */
        write_checkpoint(interval + 1, n_stop, interval_size,
                         orig_interval_start, orig_initial_size,
                         orig_intervals, effective_intervals, max_pidx);
    }

    if (csv) fclose(csv);

    /* Remove checkpoint on clean completion — the run is done */
    remove(CHECKPOINT_FILE);

    printf("\nDone.\n");
    return 0;
}

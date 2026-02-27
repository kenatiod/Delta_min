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

    int table_size = max_pidx + 50; 
    if (table_size > MAX_PRIMES) table_size = MAX_PRIMES;
    build_prime_table(table_size);
    
    if (max_pidx > num_primes) max_pidx = num_primes;

    int nthreads = omp_get_max_threads();
    char buf1[64];

    printf("\nDelta_min: Formatted Output & Optimized Logic\n");
    printf("=================================================\n\n");
    printf("Tracking: min(Pidx-omega), max(omega), min(Pidx) per interval.\n");
    printf("Intervals: %d doubling, initial size %s\n", intervals, fmt_num(interval_size, buf1, sizeof(buf1)));
    printf("Threads: %d, max prime index: %d (p_%d = %lu)\n\n",
           nthreads, max_pidx, max_pidx, (unsigned long)primes[max_pidx]);

 
    printf("            Interval Start  minDelta = Pidx - omega      n@minDelta  maxOmmega     n@maxOmega  minPidx        n@minPidx     time\n");
    printf("--------------------------------------------------------------------------------------------------------------------------------------\n");

    if (interval_start == 1 && interval_size == 1) { // Special line inserted for n=1, if that is the start
        printf("                         1         0     0     0                  0      0                  1      0                  1    0.00s\n");
        interval_start = 2;
        interval_size = 2;
        intervals--;
    }
    fflush(stdout);
    uint64_t n_start, n_stop;
    FILE *csv = fopen("delta_min_data.csv", "w");
    if (csv) {
        fprintf(csv, "interval,n_start,n_stop,min_delta,md_pidx,md_omega,md_n,max_omega,mo_n,min_pidx,mp_n\n");
    }

    for (int interval = 0; interval < intervals; interval++) {
        if (interval == 0) {
            n_start = interval_start;
        } else {
            n_start = n_stop;
        }
        n_stop = n_start + interval_size;
        interval_size *= 2;

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
            if (!find_pidx_omega(n_start, g_min_delta + 22, &p1, &o1)){
                p1 = 9999999;
                o1 = 0;
            }
            #pragma omp for schedule(dynamic, 65536)           
            for (uint64_t n = n_start+1; n <= n_stop; n++) {
                
                // Dynamic Pruning
                int current_best = g_min_delta; 
                int limit = current_best + 22; 
                if (limit > max_pidx) limit = max_pidx;

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
    }

    if (csv) fclose(csv);
    printf("\nDone.\n");
    return 0;
}

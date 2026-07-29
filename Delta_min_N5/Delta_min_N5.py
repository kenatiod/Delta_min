# Delta_min_N5.py
"""
This is a Python implementation of the Delta_Nx_min.c family of programs for N_5 = n(n+1)(n+2)(n+3)(n+4).
The parameters:
        intervals = 20
        interval_start = 1
        interval_size = 1
are written into the code because a single run is going to tell the story. Using Python's big number
capabilities, the complete product N_5 is computed for every n, and then sympy provides the
prime factor list, from which, Iota and omega are calculated. The point of this progrem
is to show that for N_5 the minIota grows so fast it outruns maxOmega right away.

By Ken Clements, May 6, 2026
"""

from sympy import prime, primefactors, primepi


intervals = 20
interval_start = 1
interval_size = 1


print(f"\nDelta_min_N5: Calculating the smallest difference between the prime number index of the greatest prime divisor\n")
if (intervals == 1):
    print(f"     of N_5 = n(n+1)(n+2)(n+3)(n+4) and the count of prime divisors of that consecutive integer product, over an interval.\n")
else:
    print(f"     of N_5 = n(n+1)(n+2)(n+3)(n+4) and the count of prime divisors of that consecutive integer product, over doubling intervals.\n")

print(f"======================================================================================================================================\n\n")

print(f"            Interval Start  minDelta = Iota - Omega      n@minDelta  maxOmega      n@maxOmega  minIota        n@minIota \n")
print(f"--------------------------------------------------------------------------------------------------------------------------------------")
if interval_start == 1:
    print(f"                         1         0     3     3                  1      3                  1      3                  1  ")
    interval_start = 2
    interval_size = 2
    intervals -= 1

for i in range(intervals):
    minDelta, Iota_at_minDelta, omega_at_minDelta, n_at_minDelta, maxOmega, n_at_maxOmega, minIota, n_at_minIota = 9999999, 0, 0, 0, 0, 0, 9999999, 0    
    for n in range(interval_start, interval_start+interval_size):
        N_5 = n * (n+1) * (n+2) * (n+3) * (n+4) 
        pf_N_5 = primefactors(N_5)
        omega = len(pf_N_5)
        Iota = int(primepi(pf_N_5[-1]))
        delta = Iota - omega
        if delta <= minDelta:
            minDelta, Iota_at_minDelta, omega_at_minDelta, n_at_minDelta = delta, Iota, omega, n
        if omega >= maxOmega:
            maxOmega, n_at_maxOmega = omega, n
        if Iota <= minIota:
            minIota, n_at_minIota = Iota, n

    print(f"{interval_start:26,} {minDelta:9,} {Iota_at_minDelta:5,} {omega_at_minDelta:5,} {n_at_minDelta:18,} {maxOmega:6,} {n_at_maxOmega:18,} {minIota:6,} {n_at_minIota:18,}")
    interval_start = interval_start+interval_size
    interval_size *= 2
print(f"\nEnd of Program.\n")



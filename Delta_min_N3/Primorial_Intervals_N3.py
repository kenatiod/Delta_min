# Primorial_Intervals_N3.py
"""
Simple program to plot the progress of iota - omega of N3 over primorial intervals,
where N3 = n(n+1)(n+2).  Starts with P_2 <= N3 < P_3 which will hold 
(P_2 = 6) <= 2*3*4=24 < (P_3 = 30). To make N fit in the intervals,
the variable, n < n_limits, runs through the cube roots of the primorials.

Uses a partial factorization over the first 50 prime numbers. Because the 
program is looking for the minimum iota, and that will be less than
50 for the ranges covered, the omega and iota values are only needed when
iota is 50 or less.

By Ken Clements 7-28-2026
"""

from sympy import primorial, prime, primerange
from math import log, exp

START_P_INTERVAL = 2
END_P_INTERVAL = 26
MAX_PRIME_INDEX = 50

IOTA_OUT_OF_RANGE = 9999

primes = [0] + list(primerange(2, prime(MAX_PRIME_INDEX+1)))

def o_and_i(n):
    omega = 0
    for i in range(1, MAX_PRIME_INDEX+1):
        iota = i
        p = primes[i]
        if n % p == 0:
            omega += 1
            while n % p == 0: n //= p
            if n == 1: return omega, iota
    return 0, IOTA_OUT_OF_RANGE
        
n, np1 = 2, 3
o_n, i_n = o_and_i(n) 
o_np1, i_np1 = o_and_i(np1) 

print(f"\n\nFor Products, N3 = n(n+1)(n+2), in Primorial Intervals, P_r <= N3 < P_(r+1):")
print(f"---+-----------------+-----------------------+----------------------+----------------------+")
print(f" r | n starts @      | Min_Delta =           | Max_Omega =          | Min_Iota =           |")
print(f"---+-----------------+-----------------------+----------------------+----------------------+")


for r in range(START_P_INTERVAL, END_P_INTERVAL+1):
    P_rp1 = primorial(r+1)
    n_start = n
    n_limit = int(exp(log(P_rp1)/3.0))            # This runs n through the cube roots of P_r
    min_d = MAX_PRIME_INDEX
    max_o = 0
    min_i = MAX_PRIME_INDEX + 1
    while n < n_limit:
        np2 = np1 + 1
        o_np2, i_np2 = o_and_i(np2)
        omega = o_n + o_np1 + o_np2
        if n % 2 == 0: omega -= 1               # Subtract 1 if even for shared factor of 2
        iota = max(i_n, i_np1, i_np2)
        if iota <= MAX_PRIME_INDEX:
            if iota <= min_i:
                min_i = iota
                n_at_min_iota = n
            if omega >= max_o:
                max_o = omega
                n_at_max_omega = n
            delta = iota - omega
            if delta <= min_d:
                min_d =  delta
                n_at_min_delta = n
        n, np1 = np1, np2
        o_n, o_np1 = o_np1, o_np2
        i_n, i_np1 = i_np1, i_np2

    print(f"{r:2} | {n_start:15,} |  {min_d:3,} @{n_at_min_delta:15,}" +
          f" | {max_o:3,} @{n_at_max_omega:15,} | {min_i:3,} @{n_at_min_iota:15,} |")

print(f"\nEnd of Program\n") 


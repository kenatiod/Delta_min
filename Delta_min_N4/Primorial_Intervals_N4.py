# Primorial_Intervals_N4.py
"""
Simple program to plot the progress of iota - omega of N4 over primorial intervals,
where N4 = n(n+1)(n+2)(n+3).  Starts with P_3 <= N4 < P_4 which will hold 
(P_3 = 30) <= 2*3*4*5=120 < (P_4 = 210). To make N4 fit in the intervals,
the variable, n < n_limits, runs through the fourth roots of the primorials.

Uses a partial factorization over the first 100 prime numbers.Because the 
program is looking for the minimum iota, and that will be less than
100 for the ranges covered, thus the omega and iota values are only needed when
iota is 100 or less.

By Ken Clements 7-28-2026
"""

from sympy import primorial, prime, primerange
from math import log, exp

START_P_INTERVAL = 3
END_P_INTERVAL = 27
MAX_PRIME_INDEX = 100

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
        
n, np1, np2 = 2, 3, 4
o_n, i_n = o_and_i(n) 
o_np1, i_np1 = o_and_i(np1) 
o_np2, i_np2 = o_and_i(np2)

print(f"\n\nFor Products, N4 = n(n+1)(n+2)(n+3), in Primorial Intervals, P_r <= N4 < P_(r+1):")
print(f"---+-----------------+-----------------------+----------------------+----------------------+")
print(f" r | n starts @      | Min_Delta =           | Max_Omega =          | Min_Iota =           |")
print(f"---+-----------------+-----------------------+----------------------+----------------------+")


for r in range(START_P_INTERVAL, END_P_INTERVAL+1):
    P_rp1 = primorial(r+1)
    n_start = n
    n_limit = round(exp(log(P_rp1)/4.0)) - 1      # This runs n through the fourth roots of P_r
    min_d = MAX_PRIME_INDEX
    max_o = 0
    min_i = MAX_PRIME_INDEX + 1
    while n < n_limit:
        np3 = np2 + 1
        o_np3, i_np3 = o_and_i(np3)
        omega = o_n + o_np1 + o_np2 + o_np3 - 1 # Subtract for shared factor of 2
        if n % 3 == 0: omega -= 1               # Subtract another on if shared factor of 3
        iota = max(i_n, i_np1, i_np2, i_np3)
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
        n, np1, np2 = np1, np2, np3
        o_n, o_np1, o_np2 = o_np1, o_np2, o_np3
        i_n, i_np1, i_np2 = i_np1, i_np2, i_np3

    print(f"{r:2} | {n_start:15,} |  {min_d:3,} @{n_at_min_delta:15,}" +
          f" | {max_o:3,} @{n_at_max_omega:15,} | {min_i:3,} @{n_at_min_iota:15,} |")

print(f"\nEnd of Program\n") 


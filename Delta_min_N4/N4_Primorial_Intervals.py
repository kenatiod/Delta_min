# M_4_Delta.py
"""
Simple program to plot the progress of iota - omega of N4 over primorial intervals,
where N4 = n(n+1)(n+2)(n+3).  Starts with P_5 <= N4 < P_6 which will hold 
(P_5 = 2310) <= 11*12*13*14=24024 < (P_6 = 30030). To make N4 fit in the intervals,
the variable, n < n_limits, runs through the fourth roots of the primorials.

Uses a partial factorization over the first 100 prime numbers.Because the 
program is looking for the minimum iota, and that will be less than
100 for the ranges covered, thus the omega and iota values are only needed when
iota is 100 or less.

By Ken Clements 7-28-2026
"""

from sympy import primorial
from math import log, exp

START_P_INTERVAL = 6
END_P_INTERVAL = 27
MAX_PRIME_INDEX = 100

primes = [0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 
    113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 
    173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 
    281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 
    349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 
    409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 
    463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 
    541]

IOTA_OUT_OF_RANGE = 9999999

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

print(f"\n\nFor Products, n_4 = n(n+1)(n+2)(n+3), in Primorial Intervals, P_r:")
for r in range(START_P_INTERVAL, END_P_INTERVAL+1):
    P_r = primorial(r)
    n_limit = round(exp(log(P_r)/4.0)) - 1      # This runs n through the fourth roots of P_r
    min_d, d_n = IOTA_OUT_OF_RANGE, IOTA_OUT_OF_RANGE
    max_o, o_n = 0, 0
    min_i, i_n = IOTA_OUT_OF_RANGE, IOTA_OUT_OF_RANGE
    while n < n_limit:
        np3 = np2 + 1
        o_np3, i_np3 = o_and_i(np3)
        omega = o_n + o_np1 + o_np2 + o_np3 - 1 # Subtract for shared factor of 2
        if n % 3 == 0: omega -= 1               # Subtract another on if shared factor of 3
        iota = max(i_n, i_np1, i_np2, i_np3)
        if iota < 100:
            if iota <= min_i:
                min_i = iota
                i_n = n
            if omega >= max_o:
                max_o = omega
                o_n = n
            delta = iota - omega
            if delta <= min_d:
                min_d =  delta
                d_n = n
        n, np1, np2 = np1, np2, np3
        o_n, o_np1, o_np2 = o_np1, o_np2, o_np3
        i_n, i_np1, i_np2 = i_np1, i_np2, i_np3

    print(f"{r=:2}  n < {n_limit:15,}  Min_Delta = {min_d:3,} @{d_n:15,}" +
          f"  Max_Omega = {max_o:3,} @{o_n:15,}   Min_Iota = {min_i:3,} @{i_n:15,}")

print(f"\nEnd of Program\n") 


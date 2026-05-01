# verify_delta_min.py
"""
Delta_min Output Verifier
=========================
Reads Delta_min_output.txt and independently checks two values per row:

  1. maxOmega  = omega(n) + omega(n+1)
                 (sum of distinct prime divisor counts of n and n+1)

  2. minPidx   = max(prime_index(gpf(n)), prime_index(gpf(n+1)))
                 (the greater prime index of the greatest prime factor
                  of n vs. n+1 — i.e. the larger of the two)

Uses sympy for reliable prime factorization of large numbers.

Usage:
    python3 verify_delta_min.py [path_to_output_file]
    (defaults to Delta_min_output.txt in the current directory)

The program also provides a good base for adding analysis routines
as it builds an internal database from the input text file.
    
By Ken Clements and Claude Sonnet 4.6, April 30, 2026    
"""

import re
import sys
from sympy import primepi, primefactors

FILE = sys.argv[1] if len(sys.argv) > 1 else "Delta_min_output.txt"

# ─── helpers ──────────────────────────────────────────────────────────────────

def omega(n):
    """Count of distinct prime divisors of n (returns 0 for n <= 1)."""
    if n <= 1:
        return 0
    return len(primefactors(n))

def Pidx(n):
    """1-based prime index of the greatest prime factor of n (0 for n <= 1)."""
    if n <= 1:
        return 0
    return int(primepi(primefactors(n)[-1]))

def verify_max_omega(n, claimed):
    computed = omega(n) + omega(n + 1)
    return computed == claimed, computed

def verify_min_pidx(n, claimed):
    Pidx_n   = Pidx(n)
    Pidx_np1 = Pidx(n + 1)
    computed = max(Pidx_n, Pidx_np1)
    return computed == claimed, computed, Pidx_n, Pidx_np1

# ─── parser ───────────────────────────────────────────────────────────────────

def strip_commas(text):
    """Remove thousands-separator commas from a string of digits."""
    prev = ""
    while prev != text:
        prev = text
        text = re.sub(r"(\d),(\d)", r"\1\2", text)
    return text

def parse_data_rows(path):
    """
    Columns (after the dashed separator):
      interval  minDelta  Pidx@minDelta  omega@minDelta  n@minDelta
      maxOmega  n@maxOmega  minPidx  n@minPidx
    Returns a list of dicts.
    """
    rows = []
    in_table = False
    with open(path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if re.match(r"^-{20,}", stripped):
                in_table = True
                continue
            if not in_table or not stripped:
                continue
            tokens = strip_commas(stripped).split()
            if len(tokens) < 9:
                continue
            try:
                rows.append({
                    "interval":       int(tokens[0]),
                    "minDelta":       int(tokens[1]),
                    "Pidx_minDelta":  int(tokens[2]),
                    "omega_minDelta": int(tokens[3]),
                    "n_minDelta":     int(tokens[4]),
                    "maxOmega":       int(tokens[5]),
                    "n_maxOmega":     int(tokens[6]),
                    "minPidx":        int(tokens[7]),
                    "n_minPidx":      int(tokens[8]),
                })
            except (ValueError, IndexError):
                continue
    return rows

# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    print(f"Reading: {FILE}")
    rows = parse_data_rows(FILE)
    print(f"Rows parsed: {len(rows)}\n")

    hdr = (f"{'Interval':>22}  {'n@maxOmega':>22}  {'claimed':>7}  {'computed':>8}  {'OK':>3}  "
           f"{'n@minPidx':>22}  {'claimed':>7}  {'computed':>8}  {'OK':>3}")
    print(hdr)
    print("─" * len(hdr))

    omega_fails, pidx_fails = [], []

    for r in rows:
        n_mo, cl_mo = r["n_maxOmega"], r["maxOmega"]
        n_mp, cl_mp = r["n_minPidx"],  r["minPidx"]

        ok_mo, comp_mo           = verify_max_omega(n_mo, cl_mo)
        ok_mp, comp_mp, ip, ip1  = verify_min_pidx(n_mp, cl_mp)

        # Cast sympy integers to plain int for formatting
        comp_mo, comp_mp, ip, ip1 = int(comp_mo), int(comp_mp), int(ip), int(ip1)

        print(
            f"{r['interval']:>22,}  {n_mo:>22,}  {cl_mo:>7}  {comp_mo:>8}  {'✓' if ok_mo else '✗':>3}  "
            f"{n_mp:>22,}  {cl_mp:>7}  {comp_mp:>8}  {'✓' if ok_mp else '✗':>3}"
        )

        if not ok_mo:
            omega_fails.append({**r, "comp": comp_mo})
        if not ok_mp:
            pidx_fails.append({**r, "comp": comp_mp, "ip": ip, "ip1": ip1})

    print()
    total_fails = len(omega_fails) + len(pidx_fails)
    if total_fails == 0:
        print("✅  All values verified correctly — no discrepancies found.")
    else:
        print(f"⚠️  {total_fails} discrepancies found  "
              f"({len(omega_fails)} maxOmega,  {len(pidx_fails)} minPidx)\n")

        if omega_fails:
            print("── maxOmega failures ──")
            for f in omega_fails:
                print(f"  Interval {f['interval']:,}:  n = {f['n_maxOmega']:,}")
                print(f"    claimed maxOmega = {f['maxOmega']},  computed = {f['comp']}")
                print(f"    omega(n) = {omega(f['n_maxOmega'])},  omega(n+1) = {omega(f['n_maxOmega']+1)}")

        if pidx_fails:
            print("── minPidx failures ──")
            for f in pidx_fails:
                print(f"  Interval {f['interval']:,}:  n = {f['n_minPidx']:,}")
                print(f"    claimed minPidx = {f['minPidx']},  computed = {f['comp']}")
                print(f"    Pidx(n) = {f['ip']},  Pidx(n+1) = {f['ip1']}")
                print(f"    max(Pidx) = {max(f['ip'], f['ip1'])}")

if __name__ == "__main__":
    main()


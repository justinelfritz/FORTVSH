#!/usr/bin/env python3
"""Compare two directories of .dat files column-by-column within a relative tolerance.

Usage: check_regression.py <reference_dir> <output_dir>
Exit 0 if all files match; exit 1 if any difference exceeds REL_TOL.
Lines starting with '#' are treated as headers and skipped.
"""
import sys
import os
import glob
import re

REL_TOL = 1e-10
ABS_FLOOR = 1e-300

def compare_dirs(ref_dir, out_dir):
    failures = []
    ref_files = sorted(glob.glob(os.path.join(ref_dir, "*.dat")))
    if not ref_files:
        print(f"ERROR: no .dat files found in {ref_dir}")
        return ["no reference files found"]

    for ref_path in ref_files:
        name = os.path.basename(ref_path)
        out_path = os.path.join(out_dir, name)
        if not os.path.exists(out_path):
            failures.append(f"{name}: output file missing")
            continue

        with open(ref_path) as rf, open(out_path) as of:
            for lineno, (rl, ol) in enumerate(zip(rf, of), 1):
                if rl.startswith('#'):
                    continue
                # Extract all numeric tokens (handles Fortran complex "(re,im)" too)
                rvals = [float(x) for x in
                         re.findall(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?', rl)]
                ovals = [float(x) for x in
                         re.findall(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?', ol)]
                for i, (rv, ov) in enumerate(zip(rvals, ovals)):
                    scale = max(abs(rv), abs(ov), ABS_FLOOR)
                    if abs(rv - ov) / scale > REL_TOL:
                        failures.append(
                            f"{name} line {lineno} col {i+1}: "
                            f"ref={rv:.6e} out={ov:.6e} "
                            f"rel_err={abs(rv-ov)/scale:.2e}")
                        if len(failures) >= 20:
                            return failures

    return failures

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <reference_dir> <output_dir>")
        sys.exit(2)

    ref_dir, out_dir = sys.argv[1], sys.argv[2]
    failures = compare_dirs(ref_dir, out_dir)

    if failures:
        for f in failures:
            print("DIFF:", f)
        sys.exit(1)

    n = len(glob.glob(os.path.join(ref_dir, "*.dat")))
    print(f"OK: all {n} reference files match within rel_tol={REL_TOL}")

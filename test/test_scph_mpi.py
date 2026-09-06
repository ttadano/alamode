#!/usr/bin/env python
"""SCPH structural optimization with the V4 array row-distributed over MPI ranks.

Every multi-rank run must reproduce the single-rank run of the same input
(rank 0 keeps the whole control flow, the other ranks only serve the V4
contractions; the partial sums differ from the serial order by roundoff only):

1. relax_kpoint : RELAX_STR = 2, default IALGO = 0 -> whole V4 slices per rank
                  (k-point builder), 1 vs 2 ranks.
2. relax_band   : RELAX_STR = 2, IALGO = 1 -> weighted unit partition (band
                  builder), 1 vs 3 ranks.
3. gamma_only   : RELAX_STR = 2 on a Gamma-only mesh with IALGO = 0 -> fewer
                  slices than ranks, the band builder is selected
                  automatically, 1 vs 2 ranks (the single-rank reference uses
                  IALGO = 1 explicitly so that both runs build V4 with the same
                  kernel; measured: the two builders agree to roundoff only,
                  which the optimization amplifies to 5.8e-8 in the final
                  displacements, 3e-10 in V0, while the same builder at 1 and 2
                  ranks agrees to 1.6e-9 and 2e-13).
4. plain_scph   : RELAX_STR = 0 (no structural optimization; the SCP solver
                  alone uses the distributed V4), 1 vs 2 ranks.
5. relax_coord  : RELAX_STR = 1 (internal coordinates only), 1 vs 2 ranks.
6. simple_mixing: IMIX = 0 (the simple-mixing SCP solver), 1 vs 2 ranks.
7. diag_only_plain : SELF_OFFDIAG = 0 (diagonal-only contraction, whole slices
                  per rank), 1 vs 2 ranks. (SELF_OFFDIAG = 0 with RELAX_STR != 0
                  is rejected by the input parser, so that combination has no case.)
8. restart_h5   : RESTART_SCPH = 1 from the h5 state of case 1 at 2 ranks: no V4
                  is built and no rank touches the service; the run must finish
                  (a hang would trip the timeout) and reproduce the state.

Tolerances: |x - y| <= max(1e-8, 1e-9 max(|x|, |y|)), the same as the other
SCPH regression tests. Not covered here: bubble corrections after a
distributed SCPH run, and a fatal error after the workers entered service
(PHON_NS::exit aborts every rank through MPI_Abort; nothing in a normal input
fails at that point).

Uses the BaTiO3 fixtures of test_batio3.py. Skipped (exit 0 with a message)
when mpirun is not available. Every run has a timeout so that a rank left
waiting in a collective fails the test instead of hanging it.
"""

import os
import shutil
import subprocess
import sys

import numpy as np

from test_batio3 import copy_input_files

PREFIX = "cBTO222_scph"
TIMEOUT = 1800


def rewrite_input(path, relax_str=None, ialgo=None, gamma_only=False, self_offdiag=None, imix=None,
                  restart=False):
    with open(path) as f:
        lines = f.readlines()
    out = []
    for line in lines:
        s = line.strip()
        if relax_str is not None and s.startswith("RELAX_STR"):
            line = "  RELAX_STR = %d\n" % relax_str
        if self_offdiag is not None and s.startswith("SELF_OFFDIAG"):
            line = "  SELF_OFFDIAG = %d\n" % self_offdiag
        if gamma_only and (s.startswith("KMESH_INTERPOLATE") or s.startswith("KMESH_SCPH")):
            line = "  %s = 1 1 1\n" % s.split("=")[0].strip()
        out.append(line)
        if s.startswith("SELF_OFFDIAG"):
            if ialgo is not None:
                out.append("  IALGO = %d\n" % ialgo)
            if imix is not None:
                out.append("  IMIX = %d\n" % imix)
            if restart:
                out.append("  RESTART_SCPH = 1\n")
    with open(path, "w") as f:
        f.writelines(out)


def run_anphon(anphonbin, input_file, logfile, nprocs):
    cmd = [anphonbin, input_file]
    if nprocs > 1:
        cmd = ["mpirun", "-np", str(nprocs)] + cmd
    try:
        with open(logfile, "w") as f:
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        print("  timeout after %d s (%s, %d ranks): suspected MPI hang" % (TIMEOUT, input_file, nprocs))
        return 1
    if proc.returncode != 0:
        print("  anphon exited with code %d (%s, %d ranks)" % (proc.returncode, input_file, nprocs))
        return 1
    return 0


def load_numbers(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    return rows


def compare_files(path_a, path_b, abs_tol=1.0e-8, rel_tol=1.0e-9):
    """|x - y| <= max(abs_tol, rel_tol * max(|x|, |y|)) for every number; every value
    must be finite and the files must hold at least one number."""
    a, b = load_numbers(path_a), load_numbers(path_b)
    if len(a) == 0 or sum(len(r) for r in a) == 0:
        print("  %s: no numeric data" % os.path.basename(path_a))
        return 1
    if len(a) != len(b):
        print("  %s: row count differs (%d vs %d)" % (os.path.basename(path_a), len(a), len(b)))
        return 1
    worst = 0.0
    for ra, rb in zip(a, b):
        if len(ra) != len(rb):
            print("  %s: column count differs" % os.path.basename(path_a))
            return 1
        for x, y in zip(ra, rb):
            if not (np.isfinite(x) and np.isfinite(y)):
                print("  %s: non-finite value" % os.path.basename(path_a))
                return 1
            if abs(x - y) > max(abs_tol, rel_tol * max(abs(x), abs(y))):
                worst = max(worst, abs(x - y))
    if worst > 0.0:
        print("  %s: max deviation %.3e beyond the tolerance" % (os.path.basename(path_a), worst))
        return 1
    return 0


def log_contains(logfile, needle):
    with open(logfile) as f:
        return needle in f.read()


def run_case(name, project_root, anphonbin, nprocs_list, relax_str=None, ialgo=None, gamma_only=False,
             expect_note=None, ialgo_reference=None, self_offdiag=None, imix=None, restart_from=None):
    """ialgo_reference: IALGO of the single-rank reference run when it must differ from
    the multi-rank runs (the two V4 builders agree to roundoff only, which the
    structural optimization amplifies to ~1e-7 in the displacements)."""
    scph_example_dir = os.path.join(project_root, "example/BaTiO3/scph_relax")
    fc_reference_dir = os.path.join(project_root, "example/BaTiO3/anharm_IFCs/4_optimize/reference")
    base = os.path.join(project_root, "test/scph_mpi")
    workdirs = {}
    for nprocs in nprocs_list:
        wd = os.path.join(base, "%s_np%d" % (name, nprocs))
        if os.path.exists(wd):
            shutil.rmtree(wd)
        os.makedirs(wd)
        os.chdir(wd)
        if copy_input_files(wd, scph_example_dir, fc_reference_dir) != 0:
            print("  could not copy the inputs for %s" % name)
            return 1
        ialgo_run = ialgo_reference if (nprocs == nprocs_list[0] and ialgo_reference is not None) else ialgo
        rewrite_input("BTO_scph_thermo.in", relax_str=relax_str, ialgo=ialgo_run, gamma_only=gamma_only,
                      self_offdiag=self_offdiag, imix=imix, restart=restart_from is not None)
        if restart_from is not None:
            # restart from the h5 state of another case: no V4 is built and no rank may
            # touch the service (the run must finish, on every rank, within the timeout)
            src = os.path.join(base, "%s_np%d" % (restart_from, nprocs_list[0]), PREFIX + ".scph.h5")
            if not os.path.exists(src):
                print("  %s: missing restart source %s" % (name, src))
                return 1
            shutil.copy(src, PREFIX + ".scph.h5")
        if run_anphon(anphonbin, "BTO_scph_thermo.in", "run.log", nprocs) != 0:
            return 1
        if restart_from is not None and not log_contains("run.log", "RESTART_SCPH is true"):
            print("  %s: the run did not restart" % name)
            return 1
        if restart_from is None and nprocs > 1 and not log_contains("run.log", "V4 rows distributed over %d MPI processes" % nprocs):
            print("  %s: the log of the %d-rank run does not report the distributed V4" % (name, nprocs))
            return 1
        if nprocs > 1 and expect_note is not None and not log_contains("run.log", expect_note):
            print("  %s: expected '%s' in the %d-rank log" % (name, expect_note, nprocs))
            return 1
        workdirs[nprocs] = wd

    files = [PREFIX + ".scph_thermo", PREFIX + ".scph_dfc2"]
    if restart_from is not None:
        # a restart skips the loop and the dfc2 output; only the postprocess is redone
        files = [PREFIX + ".scph_thermo"]
    elif relax_str != 0:
        files += [PREFIX + ".V0", PREFIX + ".normal_disp", PREFIX + ".atom_disp"]
        if relax_str != 1:
            files += [PREFIX + ".umn_tensor"]
    ref = workdirs[nprocs_list[0]]
    info = 0
    for nprocs in nprocs_list[1:]:
        for fname in files:
            info += compare_files(os.path.join(ref, fname), os.path.join(workdirs[nprocs], fname))
    return 1 if info > 0 else 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=str, default="all", help="one case name, or all")
    args = parser.parse_args()

    build_dir = os.getcwd()
    project_root = os.path.dirname(build_dir)
    anphonbin = "%s/_build/anphon/anphon" % project_root

    if shutil.which("mpirun") is None:
        print("SCPH MPI: mpirun not found, skipped")
        sys.exit(0)

    cases = [
        ("relax_kpoint", dict(nprocs_list=[1, 2])),
        ("relax_band", dict(nprocs_list=[1, 3], ialgo=1)),
        ("gamma_only", dict(nprocs_list=[1, 2], gamma_only=True, ialgo_reference=1,
                            expect_note="the band-parallel builder (IALGO = 1) is used")),
        ("plain_scph", dict(nprocs_list=[1, 2], relax_str=0)),
        ("relax_coord", dict(nprocs_list=[1, 2], relax_str=1)),
        ("simple_mixing", dict(nprocs_list=[1, 2], imix=0)),
        ("diag_only_plain", dict(nprocs_list=[1, 2], relax_str=0, self_offdiag=0)),
        ("restart_h5", dict(nprocs_list=[1, 2], restart_from="relax_kpoint")),
    ]
    known = [name for name, _ in cases]
    if args.case != "all" and args.case not in known:
        print("unknown case %s (choose from %s)" % (args.case, ", ".join(known)))
        sys.exit(1)
    failed = 0
    for name, kwargs in cases:
        if args.case != "all" and args.case != name:
            continue
        info = run_case(name, project_root, anphonbin, **kwargs)
        print("SCPH MPI [%s] --> %s" % (name, "pass" if info == 0 else "failed"))
        failed += info

    if failed == 0:
        print("SCPH MPI --> pass")
        sys.exit(0)
    print("SCPH MPI --> failed")
    sys.exit(1)

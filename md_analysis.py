#!/usr/bin/env python3
"""
GROMACS MD trajectory analysis: RMSD, RMSF, Rg and SASA.

Workflow (same as the original script):
    1. Reduce trajectory frames (gmx trjconv -skip)
    2. Check trajectory (gmx check)
    3. Remove PBC and center (gmx trjconv -pbc mol -center)
    4. Fit trajectory (gmx trjconv -fit rot+trans)
    5. RMSD, RMSF, Rg, SASA (gmx rms / rmsf / gyrate / sasa) + plots

Pixi Usage:
    # Run full analysis with custom prefix (pass script arguments after '--'):
    pixi run analyze -- -p "oxa-232"

    # Pass additional options to Pixi tasks:
    pixi run analyze -- --skip 10 -p "oxa-232"
    pixi run rmsd -- -p "oxa-232"
    pixi run quick-test -- -n index.ndx

Examples (direct Python execution):
    # Defaults (md_0_1.xtc / md_0_1.tpr in the current directory)
    python md_analysis.py --prefix MyProtein

    # Custom input files and no frame skipping
    python md_analysis.py -f run.xtc -s run.tpr --prefix OXA-48 --skip 1

    # Only RMSD and RMSF, with a custom index file and output folders
    python md_analysis.py -f md.xtc -s md.tpr -n index.ndx --prefix mutant \
        --analyses rmsd rmsf --temp-dir tmp --figures-dir plots

    # Highlight residues 144-156 on the RMSF plot, show plots interactively
    python md_analysis.py --prefix MyProtein --analyses rmsf --highlight 144 156 --show

    # MPI build of GROMACS
    python md_analysis.py --prefix MyProtein --gmx gmx_mpi
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np

ANALYSES = ("rmsd", "rmsf", "rg", "sasa")


class GromacsError(RuntimeError):
    """Raised when a GROMACS command cannot be run or fails."""


# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run GROMACS RMSD, RMSF, Rg and SASA analysis and plot the results.",
        epilog=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    io = parser.add_argument_group("input / output")
    io.add_argument("-f", "--xtc", default="md_0_1.xtc", help="input trajectory (default: %(default)s)")
    io.add_argument("-s", "--tpr", default="md_0_1.tpr", help="run input file (default: %(default)s)")
    io.add_argument("-n", "--index", default=None, help="optional index file (.ndx) passed to GROMACS tools")
    io.add_argument("-p", "--prefix", default="MyProtein",
                    help="prefix added to all generated files (default: %(default)s)")
    io.add_argument("--temp-dir", default="temp",
                    help="directory for auxiliary files: trajectories and .xvg data (default: %(default)s)")
    io.add_argument("--figures-dir", default="figures",
                    help="directory for output figures (.png) (default: %(default)s)")

    run = parser.add_argument_group("analysis options")
    run.add_argument("--gmx", default="gmx", help="GROMACS executable, e.g. gmx_mpi (default: %(default)s)")
    run.add_argument("--analyses", nargs="+", choices=ANALYSES, default=list(ANALYSES),
                     help="analyses to run (default: all)")
    run.add_argument("--skip", type=int, default=5,
                     help="keep every Nth frame; use 1 to disable (default: %(default)s)")
    run.add_argument("--time-unit", default="ns", choices=["fs", "ps", "ns", "us", "ms", "s"],
                     help="time unit for RMSD/Rg/SASA (default: %(default)s)")

    groups = parser.add_argument_group("index group selections (numbers/names from the GROMACS prompt)")
    groups.add_argument("--system-group", default="0", help="group for output/skip step (default: %(default)s)")
    groups.add_argument("--center-group", default="3", help="group used to center the protein (default: %(default)s)")
    groups.add_argument("--fit-group", default="3", help="group used for least-squares fit (default: %(default)s)")
    groups.add_argument("--rmsd-group", default="3", help="group for RMSD calculation (default: %(default)s)")
    groups.add_argument("--rmsf-group", default="3", help="group for RMSF calculation (default: %(default)s)")
    groups.add_argument("--rg-group", default="3", help="group for radius of gyration (default: %(default)s)")
    groups.add_argument("--sasa-group", default="1", help="group for SASA (default: %(default)s)")

    plot = parser.add_argument_group("plotting")
    plot.add_argument("--show", action="store_true", help="display plots interactively (default: save only)")
    plot.add_argument("--dpi", type=int, default=300, help="figure resolution (default: %(default)s)")
    plot.add_argument("--highlight", nargs=2, type=int, metavar=("START", "END"),
                      help="highlight a residue range on the RMSF plot")
    plot.add_argument("-v", "--verbose", action="store_true", help="print GROMACS output")

    return parser.parse_args(argv)


# ----------------------------------------------------------------------------
# Validation and environment checks
# ----------------------------------------------------------------------------
def find_gromacs(executable):
    """Return the full path of the GROMACS executable or raise a clear error."""
    path = shutil.which(executable)
    if path is None:
        raise GromacsError(
            f"GROMACS executable '{executable}' was not found in your PATH.\n"
            "  - Load/install GROMACS (e.g. 'module load gromacs' or 'conda install -c conda-forge gromacs')\n"
            "  - or source GMXRC:  source /usr/local/gromacs/bin/GMXRC\n"
            "  - or pass the executable name with --gmx (e.g. --gmx gmx_mpi)"
        )
    return path


def validate_args(args):
    """Check input files and argument values before running anything."""
    missing = [p for p in (args.xtc, args.tpr, args.index) if p and not Path(p).is_file()]
    if missing:
        raise FileNotFoundError("Input file(s) not found: " + ", ".join(missing))
    if args.skip < 1:
        raise ValueError("--skip must be 1 or greater")
    if args.highlight and args.highlight[0] > args.highlight[1]:
        raise ValueError("--highlight START must be <= END")


def make_directories(*directories):
    """Create each directory (including parents) if it does not exist; return them as Paths."""
    paths = []
    for directory in directories:
        path = Path(directory)
        if path.exists():
            print(f"{path} directory already exists")
        else:
            path.mkdir(parents=True)
            print(f"{path} directory created successfully")
        paths.append(path)
    return paths


# ----------------------------------------------------------------------------
# GROMACS execution
# ----------------------------------------------------------------------------
def run_gmx(gmx, subcommand, options, selections=None, verbose=False):
    """
    Run `gmx <subcommand> <options>`, feeding interactive group selections via stdin.

    selections: list of group numbers/names, equivalent to `echo 3 0 | gmx ...`.
    """
    cmd = [gmx, subcommand] + [str(o) for o in options]
    stdin_text = "\n".join(str(s) for s in selections) + "\n" if selections else None
    print("Running:", " ".join(cmd), f"(selections: {' '.join(map(str, selections))})" if selections else "")

    try:
        result = subprocess.run(cmd, input=stdin_text, capture_output=True, text=True, check=True)
    except FileNotFoundError as exc:
        raise GromacsError(f"Could not execute '{gmx}': {exc}") from exc
    except subprocess.CalledProcessError as exc:
        raise GromacsError(
            f"'gmx {subcommand}' failed with exit code {exc.returncode}.\n"
            f"Command: {' '.join(cmd)}\n--- GROMACS output ---\n{exc.stderr.strip()}"
        ) from exc

    if verbose:
        print(result.stdout)
        print(result.stderr)
    return result


def index_option(index):
    return ["-n", index] if index else []


def skip_frames(gmx, args, out_file):
    """Reduce number of frames (original: -skip 5). Returns trajectory to use."""
    if args.skip == 1:
        print("Frame skipping disabled; using the original trajectory.")
        return args.xtc
    run_gmx(gmx, "trjconv",
            ["-f", args.xtc, "-s", args.tpr, "-skip", args.skip, "-o", out_file] + index_option(args.index),
            [args.system_group], args.verbose)
    return str(out_file)


def check_trajectory(gmx, xtc, verbose):
    """Print trajectory details (number of frames, time range, ...)."""
    result = run_gmx(gmx, "check", ["-f", xtc], verbose=False)
    print(result.stderr.strip() or result.stdout.strip())


def remove_pbc(gmx, args, xtc, out_file):
    run_gmx(gmx, "trjconv",
            ["-s", args.tpr, "-f", xtc, "-o", out_file, "-pbc", "mol", "-center"] + index_option(args.index),
            [args.center_group, args.system_group], args.verbose)


def fit_trajectory(gmx, args, xtc, out_file):
    """Rotational + translational least-squares fit (needed as input for the analyses below)."""
    run_gmx(gmx, "trjconv",
            ["-s", args.tpr, "-f", xtc, "-o", out_file, "-fit", "rot+trans"] + index_option(args.index),
            [args.fit_group, args.system_group], args.verbose)


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def load_xvg(path):
    """Load a two-column (or more) .xvg file, skipping '#' and '@' comment lines."""
    try:
        return np.loadtxt(path, comments=("#", "@"))
    except OSError as exc:
        raise FileNotFoundError(f"Expected GROMACS output not found: {path}") from exc


def print_stats(y):
    print("Mean= ", np.round(np.mean(y), 2))
    print("STD= ", np.round(np.std(y), 2))


def plot_series(x, y, label, xlabel, ylabel, out_file, dpi, show, legend=True, highlight=None):
    import matplotlib.pyplot as plt  # imported here so the backend is set first

    plt.figure(figsize=(8, 5))
    plt.plot(x, y, label=label, color="b", linewidth=1)
    if highlight:
        plt.axvspan(highlight[0], highlight[1], color="green", alpha=0.1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend:
        plt.legend()
    plt.savefig(out_file, dpi=dpi)
    print(f"Saved figure: {out_file}")
    if show:
        plt.show()
    plt.close()


# ----------------------------------------------------------------------------
# Individual analyses
# ----------------------------------------------------------------------------
def analyze_rmsd(gmx, args, fitted, temp, figs):
    xvg = temp / f"{args.prefix}_rmsd.xvg"
    # Original used "echo 3 3": group for least-squares fit, then group for RMSD
    run_gmx(gmx, "rms", ["-s", args.tpr, "-f", fitted, "-o", xvg, "-tu", args.time_unit] + index_option(args.index),
            [args.fit_group, args.rmsd_group], args.verbose)
    data = load_xvg(xvg)
    time, rmsd = data[:, 0], data[:, 1]
    plot_series(time, rmsd, f"RMSD({args.prefix})", f"Time ({args.time_unit})", "RMSD (nm)",
                figs / f"{args.prefix}_RMSD.png", args.dpi, args.show)
    print_stats(rmsd)


def analyze_rmsf(gmx, args, fitted, temp, figs):
    xvg = temp / f"{args.prefix}_rmsf.xvg"
    run_gmx(gmx, "rmsf", ["-s", args.tpr, "-f", fitted, "-o", xvg, "-res"] + index_option(args.index),
            [args.rmsf_group], args.verbose)
    data = load_xvg(xvg)
    resnum, rmsf = data[:, 0], data[:, 1]
    plot_series(resnum, rmsf, f"RMSF ({args.prefix})", "Residue Number", "RMSF (nm)",
                figs / f"{args.prefix}_RMSF.png", args.dpi, args.show, highlight=args.highlight)


def analyze_rg(gmx, args, fitted, temp, figs):
    xvg = temp / f"{args.prefix}_rg.xvg"
    run_gmx(gmx, "gyrate", ["-s", args.tpr, "-f", fitted, "-o", xvg, "-tu", args.time_unit] + index_option(args.index),
            [args.rg_group], args.verbose)
    data = load_xvg(xvg)
    time, rg = data[:, 0], data[:, 1]
    plot_series(time, rg, f"Rg ({args.prefix})", f"Time({args.time_unit})", "Rg(nm)",
                figs / f"{args.prefix}_Rg.png", args.dpi, args.show)
    print_stats(rg)


def analyze_sasa(gmx, args, fitted, temp, figs):
    xvg = temp / f"{args.prefix}_sasa.xvg"
    run_gmx(gmx, "sasa", ["-s", args.tpr, "-f", fitted, "-o", xvg, "-tu", args.time_unit] + index_option(args.index),
            [args.sasa_group], args.verbose)
    data = load_xvg(xvg)
    time, sasa = data[:, 0], data[:, 1]
    plot_series(time, sasa, f"SASA ({args.prefix})", f"Time({args.time_unit})", "SASA($nm^2$)",
                figs / f"{args.prefix}_SASA.png", args.dpi, args.show, legend=False)
    print_stats(sasa)


ANALYSIS_FUNCTIONS = {
    "rmsd": analyze_rmsd,
    "rmsf": analyze_rmsf,
    "rg": analyze_rg,
    "sasa": analyze_sasa,
}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)

    if not args.show:
        matplotlib.use("Agg")  # no display needed (e.g. on HPC nodes)

    try:
        gmx = find_gromacs(args.gmx)
        validate_args(args)
        temp, figs = make_directories(args.temp_dir, args.figures_dir)

        # Preprocess trajectory
        skipped = skip_frames(gmx, args, temp / f"{args.prefix}_skipped.xtc")
        check_trajectory(gmx, skipped, args.verbose)
        no_pbc = temp / f"{args.prefix}_noPBC.xtc"
        fitted = temp / f"{args.prefix}_fitted.xtc"
        remove_pbc(gmx, args, skipped, no_pbc)
        fit_trajectory(gmx, args, str(no_pbc), fitted)

        # Run selected analyses (in original order)
        for name in ANALYSES:
            if name in args.analyses:
                ANALYSIS_FUNCTIONS[name](gmx, args, str(fitted), temp, figs)

    except (GromacsError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        return 130

    print("Analysis completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

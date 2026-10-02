# Molecular Dynamics (MD) Trajectory Analysis with Pixi & GROMACS

This project provides a reproducible, Pixi-managed environment and automated workflow for GROMACS molecular dynamics (MD) trajectory analysis. It automates trajectory preprocessing (frame reduction, PBC removal, rotational/translational fitting) and computes key biophysical properties:
- **RMSD** (Root Mean Square Deviation)
- **RMSF** (Root Mean Square Fluctuation)
- **Rg** (Radius of Gyration)
- **SASA** (Solvent Accessible Surface Area)

It also generates publication-quality plots (`.png`) and summary statistics (mean and standard deviation).

---

## 1. Project Structure

```text
MD_analysis/
├── pixi.toml            # Pixi manifest (dependencies, platforms, tasks)
├── pixi.lock            # Fully reproducible lockfile for multi-platform environments
├── .gitignore           # Git ignore rules for Pixi envs, outputs, and GROMACS backups
├── md_analysis.py       # Main Python MD analysis and plotting script
├── index.ndx            # GROMACS index file containing atom groups
├── md_0_1.tpr           # GROMACS run input / topology file
├── md_0_1.xtc           # Input trajectory file
├── md_0_1.gro           # Coordinate structure file
├── figures/             # Generated publication figures (.png)
│   ├── MyProtein_RMSD.png
│   ├── MyProtein_RMSF.png
│   ├── MyProtein_Rg.png
│   └── MyProtein_SASA.png
├── temp/                # Preprocessed trajectories and raw GROMACS XVG data
│   ├── MyProtein_skipped.xtc
│   ├── MyProtein_noPBC.xtc
│   ├── MyProtein_fitted.xtc
│   ├── MyProtein_rmsd.xvg
│   ├── MyProtein_rmsf.xvg
│   ├── MyProtein_rg.xvg
│   └── MyProtein_sasa.xvg
└── README.md            # Project documentation
```

---

## 2. Dependencies & Environment Management

All dependencies (including the GROMACS binary suite and Python scientific stack) are fully managed via [Pixi](https://pixi.sh/) using the `conda-forge` channel. No system-level installation of GROMACS or manual compilation is needed.

### Locked Dependencies
- **GROMACS** (`>=2026.3,<2027`): Full command-line suite (`gmx`, `gmx check`, `gmx trjconv`, `gmx rms`, `gmx rmsf`, `gmx gyrate`, `gmx sasa`).
- **Python** (`>=3.14.7,<3.15`): Execution runtime for the analysis and plotting wrapper.
- **NumPy** (`>=2.5.3,<3`): Data extraction and statistical calculations (mean, standard deviation).
- **Matplotlib** (`>=3.11.2,<4`): Plotting engine (automatically runs headless with `Agg` backend when `--show` is omitted).

### Supported Platforms
- `linux-64`
- `osx-64`
- `osx-arm64` (Apple Silicon)

*(Note: GROMACS binaries on conda-forge are POSIX-only. Windows users should run this project inside WSL / WSL2).*

---

## 3. Installation & Quick Start

### Prerequisites
Install [Pixi](https://pixi.sh/) on your system:
```bash
# Linux / macOS / WSL
curl -fsSL https://pixi.sh/install.sh | bash
```

### Setup Environment
Clone or navigate into the project directory and install the environment:
```bash
cd MD_analysis
pixi install
```
Pixi resolves the exact dependencies from `pixi.lock` and creates a self-contained environment in `.pixi/envs/default`.

---

## 4. Execution & Available Tasks

Pixi tasks provide one-command shortcuts to run the entire pipeline or individual analyses:

### Run the Full Analysis Pipeline
Runs frame reduction (`skip=5`), PBC correction, fitting, all 4 analyses, and generates all plots:
```bash
pixi run analyze
```

### Quick Verification Test
Performs a fast smoke test using a larger skip factor (`--skip 50`):
```bash
pixi run quick-test
```

### Run Individual Analyses
To run specific analyses only:
```bash
pixi run rmsd        # RMSD only
pixi run rmsf        # RMSF only
pixi run rg          # Radius of Gyration only
pixi run sasa        # SASA only
```

### Check Trajectory Integrity
Inspects the trajectory frames, time step, and atom count with GROMACS:
```bash
pixi run check-traj
```

### Clean Output Files
Removes generated `temp/` and `figures/` directories:
```bash
pixi run clean
```

### Display Help and Options
Prints script usage, arguments, and Pixi execution examples:
```bash
pixi run help
```

### Interactive Shell
To inspect the environment or run GROMACS commands interactively:
```bash
pixi shell
gmx --version
exit
```

---

## 5. Custom Command-Line Options

You can pass custom arguments to predefined Pixi tasks using the standard `--` separator, or run `md_analysis.py` directly:

### Passing Arguments to Pixi Tasks (Recommended)
```bash
# Set custom protein prefix on full analysis
pixi run analyze -- -p "oxa-232"

# Adjust frame skipping (e.g. skip 10 frames) with custom prefix
pixi run analyze -- --skip 10 -p "oxa-232"

# Run specific analysis with custom prefix
pixi run rmsd -- -p "oxa-232"

# Pass a custom index file
pixi run quick-test -- -n index.ndx
```

### Direct Script Execution via Pixi
```bash
# Custom prefix, trajectory, and TPR
pixi run python md_analysis.py -f md_0_1.xtc -s md_0_1.tpr -p MyProtein

# Pass custom index file
pixi run python md_analysis.py -n index.ndx

# No frame skipping (use full trajectory)
pixi run python md_analysis.py --skip 1

# Highlight a specific residue range (e.g. active site loop 144-156) on RMSF plot
pixi run python md_analysis.py --analyses rmsf --highlight 144 156

# High-resolution output (e.g. 600 DPI)
pixi run python md_analysis.py --dpi 600

# View all options
pixi run help
# or
pixi run python md_analysis.py --help
```

---

## 6. Input Requirements

The workflow expects standard GROMACS files in the working directory (or paths specified via CLI flags):
1. **Trajectory file (`.xtc`)**: Coordinates over time (e.g., `md_0_1.xtc`).
2. **Run input / topology file (`.tpr`)**: Contains topology, atom types, and default group definitions (e.g., `md_0_1.tpr`).
3. **Index file (`.ndx`)** *(Optional)*: Custom index groups if non-standard selections are needed (e.g., `index.ndx`).
4. **Group Selections**: Default group indices follow standard GROMACS conventions:
   - Group `0`: `System` (used for output trajectory)
   - Group `1`: `Protein` (used for SASA calculation)
   - Group `3`: `C-alpha` (used for centering, fitting, RMSD, RMSF, and Rg)

---

## 7. Expected Outputs & Verification

Upon running `pixi run analyze`, the following outputs are generated:

### Figures (`figures/`)
- `MyProtein_RMSD.png`: Backbone/C-alpha RMSD over time (nm).
- `MyProtein_RMSF.png`: Per-residue C-alpha root mean square fluctuation (nm).
- `MyProtein_Rg.png`: Radius of gyration over time (nm).
- `MyProtein_SASA.png`: Total solvent accessible surface area over time (nm²).

### Data Files (`temp/`)
- Preprocessed trajectories: `*_skipped.xtc`, `*_noPBC.xtc`, `*_fitted.xtc`.
- Raw Grace data files: `*_rmsd.xvg`, `*_rmsf.xvg`, `*_rg.xvg`, `*_sasa.xvg`.

### Typical Output Statistics (Test Dataset: `MyProtein`, 500 ns trajectory, `skip=5`):
- **RMSD**: Mean ≈ `0.31 nm`, STD ≈ `0.06 nm`
- **Rg**: Mean ≈ `1.78 nm`, STD ≈ `0.02 nm`
- **SASA**: Mean ≈ `126.61 nm²`, STD ≈ `3.10 nm²`

---

## 8. Limitations & Notes

- **Platform Limitation**: Conda-forge does not distribute native Windows builds of GROMACS. Windows users must use Windows Subsystem for Linux (WSL / WSL2) to run this Pixi environment.
- **Disk Space**: Trajectory preprocessing writes intermediate trajectories (`_skipped.xtc`, `_noPBC.xtc`, `_fitted.xtc`). Ensure sufficient disk space is available (approx. 1 GB for 2000 frames of a ~46k atom system). Use `pixi run clean` to remove them when finished.
- **Display Backend**: Matplotlib uses the headless `Agg` backend by default so runs do not require an active X11/GUI display. Use `--show` if running interactively with a display server.

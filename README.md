# GROMACS Molecular Dynamics (MD) Trajectory Analysis 

This respiratory provides a  automated workflow for GROMACS molecular dynamics (MD) trajectory analysis. It automates trajectory preprocessing (frame reduction, PBC removal, rotational/translational fitting) and computes key biophysical properties:
- **RMSD** (Root Mean Square Deviation)
- **RMSF** (Root Mean Square Fluctuation)
- **Rg** (Radius of Gyration)
- **SASA** (Solvent Accessible Surface Area)

It also generates publication-quality plots (`.png`) .
## 1. Dependencies & Environment Management

All dependencies (including the GROMACS binary suite and Python scientific stack) are fully managed via [Pixi](https://pixi.sh/) using the `conda-forge` channel. No system-level installation of GROMACS or manual compilation is needed.
### Dependencies
- **GROMACS** (`>=2026.3,<2027`)
- **Python** (`>=3.14.7,<3.15`)
- **NumPy** (`>=2.5.3,<3`)
- **Matplotlib** (`>=3.11.2,<4`)

#### Supported Platforms
- `linux-64`
- `osx-64`
- `osx-arm64` (Apple Silicon)

*(Note: GROMACS binaries on conda-forge are POSIX-only. Windows users should run this project inside WSL / WSL2).*

---

## 2. Installation & Quick Start

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

## 3. Execution & Available Tasks

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

## 4. Custom Command-Line Options

You can pass custom arguments to predefined Pixi tasks using the standard `--` separator, or run `md_analysis.py` directly:

### Passing Arguments to Pixi Tasks (Recommended)
```bash
# Set custom protein prefix on full analysis
pixi run analyze -- -p "my_protein" #"my_protein" can be any name that you want to give for your system.

# Adjust frame skipping (e.g. skip 10 frames) with custom prefix
pixi run analyze -- --skip 10 -p "my_protein"

# Run specific analysis with custom prefix
pixi run rmsd -- -p "my_protein"

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

## 5. Input Requirements

The workflow expects standard GROMACS files in the working directory (or paths specified via CLI flags):
1. **Trajectory file (`.xtc`)**: Coordinates over time (e.g., `md_0_1.xtc`).
2. **Run input / topology file (`.tpr`)**: Contains topology, atom types, and default group definitions (e.g., `md_0_1.tpr`).
3. **Index file (`.ndx`)** *(Optional)*: Custom index groups if non-standard selections are needed (e.g., `index.ndx`).
4. **Group Selections**: Default group indices follow standard GROMACS conventions:
   - Group `0`: `System` (used for output trajectory)
   - Group `1`: `Protein` (used for SASA calculation)
   - Group `3`: `C-alpha` (used for centering, fitting, RMSD, RMSF, and Rg)

---

## 6. Expected Outputs & Verification

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

## 7. Limitations & Notes

- **Platform Limitation**: Conda-forge does not distribute native Windows builds of GROMACS. Windows users must use Windows Subsystem for Linux (WSL / WSL2) to run this Pixi environment.
- **Disk Space**: Trajectory preprocessing writes intermediate trajectories (`_skipped.xtc`, `_noPBC.xtc`, `_fitted.xtc`). Ensure sufficient disk space is available (approx. 1 GB for 2000 frames of a ~46k atom system). Use `pixi run clean` to remove them when finished.
- **Display Backend**: Matplotlib uses the headless `Agg` backend by default so runs do not require an active X11/GUI display. Use `--show` if running interactively with a display server.

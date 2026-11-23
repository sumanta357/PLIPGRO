# PLIPGRO

PLIPGRO is a workflow for profiling protein–ligand and protein–metal interactions across GROMACS molecular dynamics trajectories.

## Requirements

- Python 3.10 or newer
- GROMACS trajectory and structure files:
  - `md_clean.xtc`
  - `md_clean.pdb`
- [PyMOL](https://pymol.org/) with Python bindings (`pymol`)
- [PLIP](https://github.com/pharmai/plip) command-line tool
- Jupyter Notebook or JupyterLab

The Python script uses only the Python standard library except for PyMOL. The notebook uses `pandas` and `matplotlib` for analysis and plotting.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/sumanta357/PLIPGRO.git
cd PLIPGRO
```

### 2. Create and activate a virtual environment

**Linux or macOS:**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows:**

```powershell
python -m venv .venv
.venv\Scripts\activate
```

### 3. Install Python packages

The repository does not currently include a dependency file. Install the packages used by the notebook and workflow:

```bash
python -m pip install --upgrade pip
python -m pip install pandas matplotlib jupyter
python -m pip install pymol
```

The PyMOL package may not be available from PyPI on every platform. If `pip install pymol` does not work, install PyMOL through Conda:

```bash
conda install -c conda-forge pymol-open-source
```

### 4. Install PLIP

Install PLIP in the same environment:

```bash
python -m pip install plip
```

Verify the command-line tools:

```bash
python -c "import pymol"
plip --version
```

## Running the trajectory analysis

1. Place `md_clean.xtc` and `md_clean.pdb` in the repository directory.
2. From the repository directory, run:

```bash
python trajectory_interaction_analysis_SK.py
```

3. Enter the ligand residue name when prompted, for example `LIG`. If the trajectory contains a metal ion to analyze, enter its residue name (for example `ZN`, `CU1`, `FE`, or `MG`); otherwise, press Enter.

The script loads each frame with PyMOL, runs PLIP, and writes a CSV summary such as:

- `LIG_md_interactions.csv`
- `LIG_ZN_md_interactions.csv`

Intermediate PDB and XML report files are removed automatically. The final CSV reports interaction signatures and the number of frames in which each interaction occurred.

## Running the notebook

Start Jupyter from the repository directory:

```bash
jupyter notebook interaction_check.ipynb
```

Alternatively, use JupyterLab:

```bash
jupyter lab
```

Update `input_file` in the notebook if the generated CSV has a name different from `LIG_NI_md_interactions.csv`.

# CTSI - Charge Transport and Signal Induction Simulator

## Overview

CTSI is a Monte Carlo simulator developed at Stanford University by **Yi Gu** (2011-2014) as part of his PhD dissertation. It simulates charge carrier drift, diffusion, and signal induction in **CZT (Cadmium Zinc Telluride)** semiconductor radiation detectors used in medical imaging (PET/SPECT).

The simulator maps high-energy photon interactions to the corresponding induced charge signals on detector electrodes, enabling accurate modeling of detector response for system design and optimization.

## Physics Modeled

- **Charge Drift**: Electron/hole transport through non-uniform electric fields
- **Diffusion**: Gaussian cloud expansion based on Einstein relation
- **Trapping**: Exponential charge loss with carrier lifetime (mu-tau product)
- **Signal Induction**: Shockley-Ramo theorem using weighting potentials
- **Noise**: Fano factor for electron-hole pair statistics, preamplifier noise

## Simulation Pipeline

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Maxwell SV │     │    GRAY     │     │    CTSI     │     │  Analysis   │
│  (E-field)  │ ──> │  (Photons)  │ ──> │ (Transport) │ ──> │ (MATLAB/Py) │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
```

1. **Maxwell SV**: 2D finite element solver generates electric field distribution
2. **GRAY**: Monte Carlo ray tracer generates photon interaction locations/energies
3. **CTSI**: Simulates charge transport and calculates electrode signals
4. **Analysis**: Post-processing with MATLAB scripts or Python ctsi_toolkit

## Directory Structure

```
CTSI/
├── ctsi/                        # Main simulation code
│   ├── *.cpp, *.h                    # C++ source files (7 modules)
│   ├── ctsi.cfg                      # Configuration file (32 parameters)
│   ├── detectorSpec.txt              # Detector geometry and material properties
│   ├── preamplifier_custom.txt       # Preamplifier filter coefficients
│   ├── event_template.txt            # Event file format reference
│   ├── compile_all.sh                # Build script
│   ├── DETECTOR_SPEC.md              # Complete technical reference
│   ├── ARCHITECTURE.md               # Pipeline and upstream format notes
│   ├── output/                       # Simulation output (gitignored)
│   └── tools/                        # Python and MATLAB tooling
│       ├── ctsi_toolkit/             # Python CLI package
│       ├── generate_efield.py        # Laplace E-field solver
│       ├── generate_uniform_efield.py # Uniform test field generator
│       ├── generate_weighting_potential.py # 3D weighting potential generator
│       ├── plot_efield.py            # E-field visualization
│       └── matlab/                   # MATLAB analysis scripts
│           ├── interactiveOut.m      # Output parser
│           └── showTraces.m          # Trajectory + waveform visualization
└── .gitignore
```

## Quick Start

### Building

```bash
cd ctsi
bash compile_all.sh
```

Or manually:
```bash
cd ctsi
g++ -c main.cpp CTSI.cpp Detector.cpp Event.cpp MyRandom.cpp Preamplification.cpp Simulation.cpp
g++ -o ctsi.exe main.o CTSI.o Detector.o Event.o MyRandom.o Preamplification.o Simulation.o
chmod a+x ctsi.exe
```

### Running

```bash
./ctsi.exe r    # List-mode: process all events from GRAY output file
./ctsi.exe i    # Interactive: manual entry of single events (x, y, z, energy)
```

### Interactive Mode Example

```
Enter the values for the event
x position (in cm, between -W/2 and W/2) = 0
y position (in cm, between -W/2 and W/2) = 0
z position (in cm, between 0 and L) = 0.25
Energy (in MeV) = 0.511
```

Output is saved to `output/interactiveOut.txt`.

### Python Toolkit

```bash
pip install numpy matplotlib

# Run a single event and plot waveforms
python -m tools.ctsi_toolkit run-event --x 0 --y 0 --z 0.25 --energy 0.511
python -m tools.ctsi_toolkit plot-waveforms output/interactiveOut.txt

# Full event visualization
python -m tools.ctsi_toolkit plot-event output/interactiveOut.txt
```

See `ctsi/tools/README.md` for full toolkit documentation.

## Data Files Required

CTSI requires pre-generated E-field and weighting potential data files that are too large for git tracking. These must be generated or obtained separately:

| File(s) | Description | How to generate |
|---------|-------------|-----------------|
| E-field grid (`.dat`) | 2D electric field from Maxwell SV or Python solver | `python tools/generate_efield.py` |
| Anode weighting potentials | 3D phi grids for each anode | `python tools/generate_weighting_potential.py` |
| Cathode weighting potentials | 3D phi grids for each cathode | `python tools/generate_weighting_potential.py` |
| GRAY event file | Photon interaction list (for list mode only) | External GRAY simulator |

Place generated files in paths matching `ctsi.cfg` lines 1-7.

## Configuration

- **`ctsi.cfg`** — 32-line config: file paths, simulation parameters, grid spacing, coordinate transforms
- **`detectorSpec.txt`** — 19 detector parameters: geometry, bias, material properties

See `ctsi/DETECTOR_SPEC.md` for complete format specifications.

## Key References

- **Yi Gu PhD Dissertation** (Stanford, 2014) — Original CTSI theory and validation
- **Medical Physics Paper** (Stanford-Hill, 2023) — High-resolution CZT PET simulation

## Authors

- **Yi Gu** — Original CTSI development (2009-2014)

# Quick Start Guide

Get CTSI running from a fresh clone in ~10 minutes.

---

## Prerequisites

| Tool | Version | Check |
|------|---------|-------|
| **g++** | Any modern version | `g++ --version` |
| **Python** | 3.10+ | `python3 --version` |
| **NumPy** | any | `pip install numpy` |
| **SciPy** | any | `pip install scipy` |
| **Matplotlib** | any | `pip install matplotlib` |

> **OS:** macOS or Linux. The binary is named `ctsi.exe` for historical reasons but is a standard Unix executable.

---

## Step 1: Build the simulator

```bash
cd ctsi/
bash compile_all.sh
```

You should see no output (no errors). Verify with:

```bash
./ctsi.exe
```

This prints the help/usage message, confirming the build worked.

---

## Step 2: Generate input data files

CTSI needs three physics input files: an electric field grid and two weighting potential grids. These are too large to include in the repo but can be generated with the included Python scripts.

### 2a. Generate the electric field

```bash
python3 tools/generate_efield.py
```

This solves the 2D Laplace equation for a 39-anode CZT detector. Takes **1-3 minutes** depending on your machine. Produces `efield_L5mm_W40mm_600V.txt`.

### 2b. Generate weighting potentials

```bash
python3 tools/generate_weighting_potential.py
```

Solves two 3D Laplace problems (anode + cathode). Takes **~1 minute**. Produces `phi_anode.txt` and `phi_cathode.txt`.

### 2c. Point the config at these files

Edit `config/ctsi.cfg` lines 1, 3, and 4 to point at the generated files:

```
1_E_field                ./efield_L5mm_W40mm_600V.txt
...
3_Anode_weighting_potential   ./phi_anode.txt
4_Cathode_weighting_potential ./phi_cathode.txt
```

Also update the grid spacing on lines 19-24 to match the weighting potential resolution. The generator uses 0.002 cm for anodes and 0.01 cm for cathodes:

```
19_Anode_phi_x_grid_space    0.002
20_Anode_phi_y_grid_space    0.002
21_Anode_phi_z_grid_space    0.002
22_Cathode_phi_x_grid_space  0.01
23_Cathode_phi_y_grid_space  0.01
24_Cathode_phi_z_grid_space  0.01
```

---

## Step 3: Run your first simulation

```bash
./ctsi.exe i
```

This starts **interactive mode**. Enter a single 511 keV photon at mid-depth:

```
Enter the values for the event
x position (in cm, between -2 and 2) = 0
y position (in cm, between -2 and 2) = 0
z position (in cm, between 0 and 0.5) = 0.25
Energy (in MeV) = 0.511
```

Wait a few seconds for the simulation to complete, then enter `n` when asked to continue.

Output is written to `output/interactiveOut.txt`.

---

## Step 4: Plot the results

```bash
python3 -m tools.ctsi_toolkit plot-waveforms output/interactiveOut.txt
```

This shows anode and cathode induced charge waveforms. You should see:
- **Anode signals**: sharp rise as electrons reach the anode plane
- **Cathode signals**: slower, opposite-polarity response from hole drift

For a full 6-panel visualization (trajectories + signals):

```bash
python3 -m tools.ctsi_toolkit plot-event output/interactiveOut.txt
```

---

## Step 5: Try the Python automation

Instead of typing into interactive mode manually, use the toolkit:

```bash
python3 -m tools.ctsi_toolkit run-event --x 0 --y 0 --z 0.25 --energy 0.511
```

Run a depth scan across the detector:

```bash
python3 -m tools.ctsi_toolkit run-zscan --x 0 --y 0 --energy 0.511 \
    --z-start 0.05 --z-end 0.45 --z-steps 9
python3 -m tools.ctsi_toolkit plot-zscan output/
```

---

## Coordinate system note

CTSI uses cathode at z=0 and anode at z=L (0.5 cm by default). When entering positions in interactive mode, z is depth from the cathode:

| z (cm) | Location |
|--------|----------|
| 0.0 | Cathode surface |
| 0.25 | Mid-depth |
| 0.5 | Anode surface |

---

## What's next

- **`DETECTOR_SPEC.md`** — full reference for all config parameters, file formats, and output layout
- **`ARCHITECTURE.md`** — how CTSI fits into the upstream simulation pipeline (Maxwell SV, GRAY)
- **`tools/README.md`** — all Python toolkit subcommands
- **`config/detectorSpec.txt`** — edit detector geometry and material properties

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `Unable to open configuration file ctsi.cfg` | You're not in the `ctsi/` directory. Run from there. |
| `Unable to open E-field file` | Check `config/ctsi.cfg` line 1 points to your generated file |
| `python3: No module named tools` | Run from `ctsi/`, not from inside `tools/` |
| `ImportError: scipy` | `pip install scipy` (needed for E-field generation only) |
| Simulation hangs | The E-field or phi grid spacing in `ctsi.cfg` doesn't match the generated files |

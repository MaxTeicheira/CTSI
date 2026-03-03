# CTSI Tools

Optional Python and MATLAB tooling for running simulations, generating input data, and visualizing results.

## Python: ctsi_toolkit

A CLI package for running simulations and plotting results. Requires Python 3.8+ with numpy and matplotlib.

```bash
# From ctsi/
pip install numpy matplotlib

# Run a single event and plot waveforms
python -m tools.ctsi_toolkit run-event --x 0 --y 0 --z 0.25 --energy 0.511
python -m tools.ctsi_toolkit plot-waveforms output/interactiveOut.txt

# Run a depth scan and plot
python -m tools.ctsi_toolkit run-zscan --x 0 --y 0 --energy 0.511 --z-start 0.05 --z-end 0.45 --z-steps 9
python -m tools.ctsi_toolkit plot-zscan output/

# Full event visualization (trajectories, signals, summary)
python -m tools.ctsi_toolkit plot-event output/interactiveOut.txt
```

### Subcommands

| Command | Description |
|---------|-------------|
| `run-event` | Run a single event via interactive mode |
| `run-zscan` | Run events at multiple depths |
| `plot-waveforms` | Plot anode/cathode waveforms from interactiveOut.txt |
| `plot-zscan` | Multi-depth waveform comparison (jet colormap) |
| `plot-event` | 6-panel event visualization |

## Python: Standalone Generators

Scripts for generating E-field and weighting potential input files.

| Script | Description |
|--------|-------------|
| `generate_efield.py` | Laplace solver for non-uniform E-field (floating BC at gaps) |
| `generate_uniform_efield.py` | Simple uniform E-field for testing |
| `generate_weighting_potential.py` | 3D weighting potential via Jacobi iteration |
| `plot_efield.py` | E-field diagnostic visualization |

```bash
# Generate E-field (outputs to current directory)
python tools/generate_efield.py

# Generate weighting potentials
python tools/generate_weighting_potential.py

# Visualize E-field
python tools/plot_efield.py
```

## MATLAB: Analysis Scripts

Located in `tools/matlab/`. Requires MATLAB.

| Script | Description |
|--------|-------------|
| `interactiveOut.m` | Parse output/interactiveOut.txt into workspace variables |
| `showTraces.m` | 3D charge trajectory + waveform visualization (calls interactiveOut.m) |

```matlab
% From ctsi/tools/matlab/
interactiveOut   % Loads and parses the data
showTraces       % Visualizes trajectories and waveforms
```

## Future Work

- Port C/A ratio depth-sensing analysis from legacy/plot_zscan.py into ctsi_toolkit

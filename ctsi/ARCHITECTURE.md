# CTSI Architecture Notes

Preserved reference information from upstream documentation.

## Simulation Pipeline

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Maxwell SV │     │    GRAY     │     │    CTSI     │     │  Analysis   │
│  (E-field)  │ ──> │  (Photons)  │ ──> │ (Transport) │ ──> │ (MATLAB/Py) │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
```

1. **Maxwell SV** (Ansoft/Ansys) — 2D finite element solver generates the electric field distribution inside the detector geometry
2. **GRAY** (Stanford) — Monte Carlo ray-driven photon transport engine; simulates Compton scattering and photoelectric absorption to produce interaction locations and energies
3. **CTSI** — Simulates charge carrier drift, diffusion, trapping, and signal induction via the Shockley-Ramo theorem
4. **Analysis** — Post-processing with MATLAB scripts or Python ctsi_toolkit

## GRAY Output Format (upstream of CTSI)

GRAY produces a whitespace-separated file with up to 12 columns. CTSI reads the first 10:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | eventType | Interaction type: 0=decay, 1=Compton, 3=photoelectric |
| 2 | eventID | Positron/event number |
| 3 | dir | Photon color/direction: 0,1=annihilation, 3=single gamma |
| 4 | time | Interaction time (seconds) |
| 5 | energy | Energy deposited (MeV) |
| 6 | x | X position (cm) |
| 7 | y | Y position (cm) |
| 8 | z | Z position (cm) |
| 9 | one | Source number (typically 1) |
| 10 | detectorID | Detector ID |
| 11* | scatterFlag | Scatter flag (not read by CTSI) |
| 12* | materialID | Material ID (not read by CTSI) |

See `DETECTOR_SPEC.md` Section 8 for the CTSI-specific event file format details and coordinate transforms.

## Maxwell SV E-field Export

When exporting electric field data from Maxwell SV for use with CTSI:

1. Export as ASCII text with 4 columns: `xPos_mm zPos_mm Ex_V/m Ez_V/m`
2. **Add a header line manually:** `numXPoints numZPoints 1 0`
3. Grid spacing must match `ctsi.cfg` parameters (lines 16-18)
4. X positions in mm, centered at 0 (range: -W/2*10 to +W/2*10)
5. Z positions in mm (range: 0 to L*10)

See `DETECTOR_SPEC.md` Section 12 for the complete file format specification.

## Alternative E-field Generation

If Maxwell SV is unavailable, the Python scripts in `tools/` can generate E-field files:
- `tools/generate_efield.py` — Solves the 2D Laplace equation with floating BCs at electrode gaps
- `tools/generate_uniform_efield.py` — Generates a simple uniform field for testing

## Key References

- **Yi Gu PhD Dissertation** (Stanford, 2014) — Original CTSI theory and validation
- **Medical Physics Paper** (Stanford-Hill, 2023) — High-resolution CZT PET simulation using CTSI

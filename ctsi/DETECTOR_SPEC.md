# CTSI Detector Specification Reference

Comprehensive reference for driving the CTSI (Charge Transport Simulation for Imaging) CZT detector simulator and parsing its results, extracted from verified source code.

---

## 1. Overview

CTSI simulates charge transport and signal induction in a **CZT (Cadmium Zinc Telluride)** semiconductor radiation detector. It uses the **Shockley-Ramo theorem** to compute induced signals on electrodes as electron and hole charge clouds drift through the detector under an applied electric field.

The simulator:
1. Reads a gamma-ray interaction event (position + energy)
2. Generates electron-hole pairs based on ionization energy and Fano statistics
3. Splits the charge cloud into configurable sub-elements
4. Transports each element through the E-field with diffusion, trapping, and detrapping
5. Computes induced current on each electrode via weighting potentials
6. Applies preamplifier shaping and noise

---

## 2. Detector Geometry

All 19 parameters from `detectorSpec.txt`:

| # | Key | Default Value | Unit | Description |
|---|-----|---------------|------|-------------|
| 1 | `L` | 0.5 | cm | Detector thickness (cathode-to-anode distance) |
| 2 | `W` | 4 | cm | Detector width (X and Y extent) |
| 3 | `NUM_ANODES` | 39 | — | Number of anode strips |
| 4 | `NUM_CATHODES` | 8 | — | Number of cathode strips |
| 5 | `ANODE_PITCH` | 1 | mm | Center-to-center anode spacing |
| 6 | `ANODE_WIDTH` | 100 | um | Physical width of each anode strip |
| 7 | `CATHODE_PITCH` | 5 | mm | Center-to-center cathode spacing |
| 8 | `CATHODE_WIDTH` | 4900 | um | Physical width of each cathode strip |
| 9 | `BIAS` | 1200 | V/cm | Anode potential relative to cathode |
| 10 | `CZT_W_FACTOR` | 4.64 | eV | Ionization energy (W-value) |
| 11 | `DIFFUSION_E` | 25.249 | cm²/s | Electron diffusion coefficient |
| 12 | `DIFFUSION_H` | 1.262 | cm²/s | Hole diffusion coefficient |
| 13 | `MU_E` | 1000 | cm²/(V·s) | Electron mobility |
| 14 | `MU_H` | 50 | cm²/(V·s) | Hole mobility |
| 15 | `TAU_E` | 1e-5 | s | Electron trapping lifetime |
| 16 | `TAU_H` | 1e-6 | s | Hole trapping lifetime |
| 17 | `FANO` | 0.089 | — | Fano factor for charge generation statistics |
| 18 | `EPSILON_R` | 10.9 | — | Relative permittivity of CZT |
| 19 | `TEMP` | 293 | K | Absolute temperature |

**Derived parameters** (computed after import in `Detector.cpp:913-914`):
- `LAMBDA_E = MU_E * BIAS * TAU_E` — Electron mean drift length (cm). Default: 1000 × 1200 × 1e-5 = **12 cm**
- `LAMBDA_H = MU_H * BIAS * TAU_H` — Hole mean drift length (cm). Default: 50 × 1200 × 1e-6 = **0.06 cm**

Diffusion coefficients are derived from the Einstein relation: `D = μ × kT/q`. For electrons at 293 K: `1000 × (1.3806503e-23 × 293) / 1.60217646e-19 = 25.25 cm²/s`.

---

## 3. Coordinate System

*Source: `Simulation.cpp:827` — "Cathode plane is at z = 0 mm, anode plane is at z = L mm."*

| Axis | Origin | Range | Notes |
|------|--------|-------|-------|
| **Z** | Cathode (z=0) | 0 to L (0.5 cm) | Anode at z=L. Electrons drift toward anode (+z), holes toward cathode (-z) |
| **X** | Detector center | -W/2 to +W/2 (-2 to +2 cm) | Anode strips indexed along X |
| **Y** | Detector center | -W/2 to +W/2 (-2 to +2 cm) | Cathode strips indexed along Y |

**All internal units are centimeters (cm).**  Pitches and widths in `detectorSpec.txt` use mm and um respectively but are converted internally.

---

## 4. Electrode Layout

*Source: `Detector.cpp:121-136` (calcOffset)*

### Anode Offsets (along X)

The anchor point for the leftmost anode is at:
```
anchor = -(NUM_ANODES/2 - 0.5) * ANODE_PITCH   [mm]
```

Each anode center in **cm**:
```
xOffset[i] = (anchor + i * ANODE_PITCH) * 0.1   [cm]
```

For default values (39 anodes, 1 mm pitch):
- Anchor = -(39/2 - 0.5) × 1 = -19 mm
- `xOffset[0]` = -19 × 0.1 = **-1.9 cm**
- `xOffset[19]` = 0 × 0.1 = **0.0 cm** (center)
- `xOffset[38]` = 19 × 0.1 = **+1.9 cm**

### Cathode Offsets (along Y)

Same formula along Y:
```
anchor = -(NUM_CATHODES/2 - 0.5) * CATHODE_PITCH   [mm]
yOffset[i] = (anchor + i * CATHODE_PITCH) * 0.1    [cm]
```

For default values (8 cathodes, 5 mm pitch):
- Anchor = -(8/2 - 0.5) × 5 = -17.5 mm
- `yOffset[0]` = -17.5 × 0.1 = **-1.75 cm**
- `yOffset[7]` = 17.5 × 0.1 = **+1.75 cm**

### Collection Logic

*Source: `Detector.cpp:74-115` (anodeCollection), `Detector.cpp:143-183` (cathodeCollection)*

**Anode collection** — an electron charge element at final X position is collected by anode `i` if:

```
|xValue*10 - floor(xValue*10 + ANODE_PITCH*0.5)| <= ANODE_WIDTH * 0.5
```

Where `xValue` is in cm, `*10` converts to mm. The charge must land within half the anode width (in um converted to mm: 100 um = 0.1 mm → `ANODE_WIDTH*0.5` = 50 um = 0.05 mm) of the nearest anode center. The anode index is computed as:

```
anodeInd = floor((xValue*10 + ANODE_PITCH*(NUM_ANODES/2)) / ANODE_PITCH)
```

**Note:** `ANODE_WIDTH` in the collection check is compared in **um** units against a value also in **um** (the expression `xValue*10 - floor(...)` produces mm-scale residuals, and `ANODE_WIDTH` is in um, so the check effectively tests whether the charge is within `ANODE_WIDTH/2` um of the anode center converted to the same scale). In the default config, `ANODE_WIDTH = 100` um and `ANODE_PITCH = 1` mm, so the collection width is `100*0.5 = 50` in the comparison units.

**Cathode collection** follows the same logic along Y using hole positions, `CATHODE_PITCH`, and `CATHODE_WIDTH`.

---

## 5. Configuration Files

### 5.1 `config/ctsi.cfg` — Run Configuration (32 lines)

Format: `N_Label<TAB>VALUE` — one parameter per line, read sequentially by line number.

| Line | Label | Type | Description | Default |
|------|-------|------|-------------|---------|
| 1 | `1_E_field` | path | E-field data file | `./currCZT/ctsi_v4_efield_xz_0p01mm.dat` |
| 2 | `2_Weighting_potential_source` | char | `n` = numeric (from file), `a` = analytic (2D formula) | `n` |
| 3 | `3_Anode_weighting_potential` | path | Anode phi data file | `./currCZT/ctsi_v4_anode_phi_0p01cm.dat` |
| 4 | `4_Cathode_weighting_potential` | path | Cathode phi data file | `./currCZT/ctsi_v4_cathode_phi_0p01cm.dat` |
| 5 | `5_GRAY_output` | path | Event input file (GRAY format) | `./single_event.txt` |
| 6 | `6_Detector_specification` | path | Detector spec file | `./config/detectorSpec.txt` |
| 7 | `7_Preamplifier_specification` | path | Preamplifier spec file | `./config/preamplifier_custom.txt` |
| 8 | `8_List-mode_output_name` | path | Output file for list-mode results | `listOut.txt` |
| 9 | `9_Random_seed` | int | PRNG seed | `65539` |
| 10 | `10_Number_of_charge_elements` | int | Sub-elements per charge cloud | `50` |
| 11 | `11_Max_num_sim_steps` | int | Max simulation time steps | `10000` |
| 12 | `12_Small_pixel_distance` | double | Small-pixel-effect distance | `0.445` |
| 13 | `13_Anode_trigger_threshold` | double | Min induced charge to trigger anode | `0.05` |
| 14 | `14_Cathode_trigger_threshold` | double | Min induced charge to trigger cathode (negative) | `-0.04` |
| 15 | `15_FWHM_preamplifier_noise` | double | FWHM of preamplifier Gaussian noise | `4.008e-16` |
| 16 | `16_E_field_x_grid_space` | double | E-field grid spacing in X (mm) | `0.01` |
| 17 | `17_E_field_y_grid_space` | double | E-field grid spacing in Y (mm) | `1` |
| 18 | `18_E_field_z_grid_space` | double | E-field grid spacing in Z (mm) | `0.01` |
| 19 | `19_Anode_phi_x_grid_space` | double | Anode phi grid spacing in X (cm) | `0.01` |
| 20 | `20_Anode_phi_y_grid_space` | double | Anode phi grid spacing in Y (cm) | `0.01` |
| 21 | `21_Anode_phi_z_grid_space` | double | Anode phi grid spacing in Z (cm) | `0.01` |
| 22 | `22_Cathode_phi_x_grid_space` | double | Cathode phi grid spacing in X (cm) | `0.01` |
| 23 | `23_Cathode_phi_y_grid_space` | double | Cathode phi grid spacing in Y (cm) | `0.01` |
| 24 | `24_Cathode_phi_z_grid_space` | double | Cathode phi grid spacing in Z (cm) | `0.01` |
| 25 | `25_Event_x_pos_scale_factor` | double | Scale factor for event X position | `1` |
| 26 | `26_Event_y_pos_scale_factor` | double | Scale factor for event Y position | `1` |
| 27 | `27_Event_z_pos_scale_factor` | double | Scale factor for event Z position | `-1` |
| 28 | `28_Event_x_pos_offset` | double | Offset added to event X position (cm) | `0` |
| 29 | `29_Event_y_pos_offset` | double | Offset added to event Y position (cm) | `0` |
| 30 | `30_Event_z_pos_offset` | double | Offset added to event Z position (cm) | `0.5` |
| 31 | `31_Output_mode` | char | `r` = pseudo-RENA format, `m` = matrix format | `r` |
| 32 | `32_Caution` | bool | Enable data consistency checks (`true`/`false`) | `false` |

### 5.2 `config/detectorSpec.txt` — Detector Parameters (19 lines)

Format: `KEY VALUE<TAB><TAB>// comment` — one parameter per line, comments are **required**.

See [Section 2](#2-detector-geometry) for the full table.

---

## 6. Z-Coordinate Transform

*Source: `main.cpp:128-132`, `configParamsStruct.h:32-37`, `CTSI.cpp:653-658`*

Positions from the event file (or interactive mode input) are transformed before simulation:

```
pos_internal = pos_input * SCALE + OFFSET
```

Applied per-axis:
```
x_internal = x_input * EVENT_POS_SCALE_X + EVENT_POS_OFFSET_X
y_internal = y_input * EVENT_POS_SCALE_Y + EVENT_POS_OFFSET_Y
z_internal = z_input * EVENT_POS_SCALE_Z + EVENT_POS_OFFSET_Z
```

**Note:** Scaling is applied **before** the offset (`configParamsStruct.h:31`).

### Default Z-transform

With default config values (`SCALE_Z = -1`, `OFFSET_Z = 0.5`):
```
z_internal = z_input * (-1) + 0.5
```

### Why this exists

The GRAY Monte Carlo output uses a coordinate system where the **anode** surface is at z=0 and depth increases away from the anode. CTSI's internal coordinate system has z=0 at the **cathode**. The default transform flips and offsets to reconcile these conventions.

### Inverse transform (for creating event files)

To convert from internal (simulator) coordinates to event-file coordinates:
```
z_event = (z_internal - OFFSET_Z) / SCALE_Z
```

With defaults:
```
z_event = (z_internal - 0.5) / (-1)
```

### Examples

| z_input (event file) | z_internal (simulator) | Physical location |
|----------------------|------------------------|-------------------|
| 0.0 | 0.5 cm | At the anode (z=L) |
| 0.25 | 0.25 cm | Mid-depth |
| 0.5 | 0.0 cm | At the cathode (z=0) |

**Interactive mode note:** In interactive mode (`main.cpp:96-149`), the user enters positions in **internal coordinates** (cm, with z between 0 and L). The transform is then applied on top of those values, so with default config the same transform applies. The prompt says `"z position (in cm, between 0 and L)"` and validates `z >= 0 && z <= 0.5`.

---

## 7. Running the Simulator

### Build

```bash
# From ctsi/
g++ -c src/Detector.cpp
g++ -c src/CTSI.cpp
g++ -c src/main.cpp
g++ -c src/Event.cpp
g++ -c src/MyRandom.cpp
g++ -c src/Preamplification.cpp
g++ -c src/Simulation.cpp
g++ -o ctsi.exe main.o CTSI.o Detector.o Event.o MyRandom.o Preamplification.o Simulation.o
chmod a+x ctsi.exe
```

Or simply: `bash compile_all.sh`

### Run Modes

| Mode | Command | Description |
|------|---------|-------------|
| **List mode** | `./ctsi.exe r` | Reads events from file specified in `ctsi.cfg` line 5, processes all events, writes to `listOut.txt` |
| **Interactive** | `./ctsi.exe i` | Prompts user for single events via stdin, writes detailed output to `interactiveOut.txt` |
| **Help** | `./ctsi.exe` | Prints help (no arguments) |

### Interactive Mode Protocol

*Source: `main.cpp:96-149`*

The simulator enters a loop prompting for one event at a time:

```
Enter the values for the event
x position (in cm, between -2 and 2) = <USER_INPUT>
y position (in cm, between -2 and 2) = <USER_INPUT>
z position (in cm, between 0 and 0.5) = <USER_INPUT>
Energy (in MeV) = <USER_INPUT>
```

Input validation:
- `x`: must be in `[-W/2, +W/2]`
- `y`: must be in `[-W/2, +W/2]`
- `z`: must be in `[0, 0.5]` (hardcoded, not `[0, L]`)
- `Energy`: must be > 0

After simulation completes, it prompts:
```
Continue? (y/n): <USER_INPUT>
```

Enter `y` for another event, `n` to quit.

**Note on energy:** The input is in MeV; it is internally converted to eV by multiplying by 1e6 (`main.cpp:129`).

**Note on transform:** Interactive mode applies the same `SCALE + OFFSET` transform to the user-entered positions (`main.cpp:130-132`). With default config, the z-transform means `z_input=0.25` → `z_internal=0.25` (midpoint).

---

## 8. Event File Format

*Source: `CTSI.cpp:590-671` (importEvents)*

The event file (GRAY output) is whitespace-separated with 10 fields per line:

```
eventType  eventID  dir  time  energy  x  y  z  one  detectorID
```

| Column | Field | Type | Description |
|--------|-------|------|-------------|
| 1 | `eventType` | int | Event type identifier |
| 2 | `eventID` | int | Event number (groups interactions) |
| 3 | `dir` | int | Photon direction index |
| 4 | `time` | float | Interaction time |
| 5 | `energy` | float | Deposited energy in **MeV** (converted to eV internally via ×1e6) |
| 6 | `x` | float | X position — in **pre-transform** coordinates |
| 7 | `y` | float | Y position — in **pre-transform** coordinates |
| 8 | `z` | float | Z position — in **pre-transform** coordinates |
| 9 | `one` | int | (Unused / always 1) |
| 10 | `detectorID` | int | Detector ID |

Events are grouped by `(eventID, dir)` pairs. When either changes, a new Event object is created. All lines with the same `(eventID, dir)` are treated as multiple interactions within a single event.

**Position transform:** All positions have `SCALE * pos + OFFSET` applied before simulation (see [Section 6](#6-z-coordinate-transform)).

---

## 9. Output: `interactiveOut.txt` Binary Layout

*Source: `CTSI.cpp:814-1077` (output function, interactive dump), cross-validated with `interactiveOut.m`*

This file is produced in interactive mode only. It is a **text file** with one number per line (not truly binary, but a sequential flat dump). All values are newline-separated ASCII doubles/ints.

### Sequential Layout

The notation uses these variables read from the data itself:
- `timeLen` — number of simulation time steps
- `numAnCollectedCharge` — number of anodes that collected charge
- `numCaCollectedCharge` — number of cathodes that collected charge
- `numAnTrig` — number of triggered anodes (met threshold)
- `numCaTrig` — number of triggered cathodes
- `numIntrxn` — number of interactions in the event
- `numChargeElems` — number of charge sub-elements per interaction
- `preampLen` — number of preamplifier time samples

```
SECTION 1: Time Vector
  timeLen                              [1 value]
  timeVec[0..timeLen-1]                [timeLen values]

SECTION 2: Collector IDs
  numAnCollectedCharge                 [1 value]
  anCollectorID[0..N-1]                [numAnCollectedCharge values: anode indices]
  numCaCollectedCharge                 [1 value]
  caCollectorID[0..N-1]                [numCaCollectedCharge values: cathode indices]

SECTION 3: Anode Trigger Data
  numAnTrig                            [1 value]
  anodeID[0..numAnTrig-1]              [numAnTrig values: triggered anode indices]
  anEnergy[0..numAnTrig-1]             [numAnTrig values: induced charge amplitude]
  anTriggerTime[0..numAnTrig-1]        [numAnTrig values: trigger timestamps]
  noisyAnEnergy[0..numAnTrig-1]        [numAnTrig values: amplitude with noise]
  noisyAnTriggerTime[0..numAnTrig-1]   [numAnTrig values: trigger time with noise]

SECTION 4: Anode Waveforms
  numAnCollectedCharge                 [1 value — count of anVsTime waveforms]
  for each anode:
    anVsTime[i][0..timeLen-1]          [timeLen values: induced current vs time]

SECTION 5: Cathode Trigger Data
  numCaTrig                            [1 value]
  cathodeID[0..numCaTrig-1]            [numCaTrig values]
  caEnergy[0..numCaTrig-1]             [numCaTrig values]
  caTriggerTime[0..numCaTrig-1]        [numCaTrig values]
  noisyCaEnergy[0..numCaTrig-1]        [numCaTrig values]
  noisyCaTriggerTime[0..numCaTrig-1]   [numCaTrig values]

SECTION 6: Cathode Waveforms
  numCaCollectedCharge                 [1 value — count of caVsTime waveforms]
  for each cathode:
    caVsTime[i][0..timeLen-1]          [timeLen values]

SECTION 7: Event Trails (charge element trajectories)
  numFieldsInTrail                     [1 value: always 10]
  numIntrxn                            [1 value]
  for each interaction i:
    numChargeElems                     [1 value]
    for each charge element j:
      len; xPosE[0..len-1]             [electron X positions]
      len; yPosE[0..len-1]             [electron Y positions]
      len; zPosE[0..len-1]             [electron Z positions]
      len; qETrapped[0..len-1]         [trapped electron charge per step]
      len; qEMobile[0..len-1]          [mobile electron charge per step]
      len; xPosH[0..len-1]             [hole X positions]
      len; yPosH[0..len-1]             [hole Y positions]
      len; zPosH[0..len-1]             [hole Z positions]
      len; qHTrapped[0..len-1]         [trapped hole charge per step]
      len; qHMobile[0..len-1]          [mobile hole charge per step]

SECTION 8: Anode Weighting Potential Traces
  numAnodes                            [1 value: count of collecting anodes]
  for each anode i:
    numIntrxn                          [1 value]
    for each interaction j:
      numChargeElems                   [1 value]
      for each charge element h:
        len; EPhi[0..len-1]            [electron weighting potential trace]
        len; HPhi[0..len-1]            [hole weighting potential trace]

SECTION 9: Cathode Weighting Potential Traces
  numCathodes                          [1 value]
  for each cathode i:
    numIntrxn                          [1 value]
    for each interaction j:
      numChargeElems                   [1 value]
      for each charge element h:
        len; EPhi[0..len-1]
        len; HPhi[0..len-1]

SECTION 10: Trail Sizes (electron)
  numIntrxn                            [1 value]
  for each interaction:
    numChargeElem                      [1 value]
    trailSizeE[i][j]                   [1 value per charge element]

SECTION 11: Trail Sizes (hole)
  numIntrxn                            [1 value]
  for each interaction:
    numChargeElem                      [1 value]
    trailSizeH[i][j]                   [1 value per charge element]

SECTION 12: Preamplifier Time Vector
  preampLen                            [1 value]
  timeVecPreamp[0..preampLen-1]        [preampLen values]

SECTION 13: Regularized Anode Signals
  for each collecting anode (numAnCollectedCharge):
    anVsTimeReg[i][0..preampLen-1]     [preampLen values]

SECTION 14: Regularized Cathode Signals
  for each collecting cathode (numCaCollectedCharge):
    caVsTimeReg[i][0..preampLen-1]     [preampLen values]

SECTION 15: Noisy Regularized Anode Signals
  for each collecting anode:
    noisyAnVsTimeReg[i][0..preampLen-1]

SECTION 16: Noisy Regularized Cathode Signals
  for each collecting cathode:
    noisyCaVsTimeReg[i][0..preampLen-1]

SECTION 17: Preamplifier Anode Signals
  for each collecting anode:
    anVsTimePreamp[i][0..preampLen-1]

SECTION 18: Preamplifier Cathode Signals
  for each collecting cathode:
    caVsTimePreamp[i][0..preampLen-1]

SECTION 19: Noisy Preamplifier Anode Signals
  for each collecting anode:
    noisyAnVsTimePreamp[i][0..preampLen-1]

SECTION 20: Noisy Preamplifier Cathode Signals
  for each collecting cathode:
    noisyCaVsTimePreamp[i][0..preampLen-1]
```

**Parsing reference:** `interactiveOut.m` (MATLAB) provides a working parser that reads this format sequentially using a running `baseInd` pointer.

---

## 10. Output: `listOut.txt` Formats

*Source: `CTSI.cpp:730-802`*

### RENA format (output mode `r`)

Header line: `Event	Channel	Charge Induced	Timestamp	Type`

Data lines (tab-separated):
```
eventCount	electrodeID	chargeInduced	triggerTime	type
```

- `type` = `A` for anode, `C` for cathode
- Cathode `chargeInduced` is negated (written as `-caEnergy[i]`)
- One line per triggered electrode per event

### Matrix format (output mode `m`)

No header line. Columns (tab-separated):

| Column | Description |
|--------|-------------|
| 1 | Event ID (sequential, incremented per event) |
| 2 | Original event ID (from input file) |
| 3 | Electrode type: `0` = anode, `1` = cathode |
| 4 | Electrode ID |
| 5 | Pulse amplitude |
| 6 | Trigger time |
| 7 | Pulse amplitude with noise |
| 8 | Trigger time with noise |

---

## 11. Preamplifier Specification

*Source: `preamplifier_custom.txt`*

Format: `N_Label<TAB><TAB>VALUE(s)` — 5 lines total.

```
1_Sampling_interval         1e-9
2_Num_anode_coefficients    2
3_Anode_coefficients        1400  -0.99999285714285714
4_Num_cathode_coefficients  2
5_Cathode_coefficients      1400  -0.99999285714285714
```

| Line | Field | Description |
|------|-------|-------------|
| 1 | Sampling interval | Time step for preamplifier output in seconds (1 ns) |
| 2 | Num anode coefficients | Number of IIR filter coefficients for anode |
| 3 | Anode coefficients | Space-separated IIR filter coefficients |
| 4 | Num cathode coefficients | Number of IIR filter coefficients for cathode |
| 5 | Cathode coefficients | Space-separated IIR filter coefficients |

The preamplifier applies an IIR (Infinite Impulse Response) filter to the raw induced current waveform. The default coefficients `[1400, -0.99999285714285714]` model a charge-sensitive preamplifier with gain and decay.

---

## 12. E-field Data File Format

*Source: `Detector.cpp:974-1076` (importEField), file header of `ctsi_v4_efield_xz_0p01mm.dat`*

The E-field file is a 2D field (X-Z plane, assumes Y-uniformity).

### Header

First line: `numXPos numZPos temp temp` — four integers.

Example: `4001 501 1 0` — 4001 X grid points, 501 Z grid points.

### Data

Remaining lines: four columns per row, whitespace-separated:

```
xPos_mm  zPos_mm  Ex_V/m  Ez_V/m
```

| Column | Description |
|--------|-------------|
| 1 | X position in mm |
| 2 | Z position in mm |
| 3 | E-field X component in V/m |
| 4 | E-field Z component in V/m |

Data is ordered with Z varying fastest (all Z values for a given X, then next X). E-field values are converted from V/m to V/cm internally by multiplying by 0.01 (`Detector.cpp:232`).

Grid spacing must match the `E_FIELD_GRID_SPACE_X` and `E_FIELD_GRID_SPACE_Z` config parameters (in mm).

---

## 13. Weighting Potential File Format

*Source: `Detector.cpp:1078-1231` (importPhi)*

Both anode and cathode weighting potentials use the same 3D format.

### Header

First line: `numX numY numZ` — three integers.

Example: `401 401 51` — 401 X points, 401 Y points, 51 Z points.

### Grid Vectors

Following the header, the grid coordinates appear sequentially:
1. **X grid** — `numX` values (in cm, exploiting symmetry so typically only half-space stored)
2. **Y grid** — `numY` values (in cm)
3. **Z grid** — `numZ` values (in cm)

### Scalar Data

After the grid vectors, the 3D scalar array is stored in **X-major, Y-middle, Z-minor** order:

```
for xIndex in 0..numX-1:
  for yIndex in 0..numY-1:
    for zIndex in 0..numZ-1:
      scalar[xIndex][yIndex][zIndex]
```

After import, all values are **normalized** so the maximum absolute value equals 1.0 (`Detector.cpp:1203-1213`).

The phi data exploits electrode symmetry — only one quadrant of the full potential is stored. During lookup, `abs(xValue)` and `abs(yValue)` are used (`Detector.cpp:399-400`), and trilinear interpolation is applied.

Grid spacing must match the phi grid space config parameters (lines 19-24 of `ctsi.cfg`), in cm.

---

## 14. Source Files & Build

### C++ Source Files

| File | Description |
|------|-------------|
| `main.cpp` | Entry point, argument parsing, interactive mode loop |
| `CTSI.cpp` / `CTSI.h` | Top-level application: config import, event import, output, simulation orchestration |
| `Detector.cpp` / `Detector.h` | Detector geometry, E-field import/interpolation, phi import/interpolation, electrode collection logic |
| `Simulation.cpp` / `Simulation.h` | Charge transport engine, signal computation (Shockley-Ramo), time stepping |
| `Event.cpp` / `Event.h` | Event data structure (interactions with position, energy, time) |
| `MyRandom.cpp` / `MyRandom.h` | Random number generator wrapper |
| `Preamplification.cpp` / `Preamplification.h` | IIR preamplifier filter, signal regularization, noise addition |

### Header-Only Files

| File | Description |
|------|-------------|
| `configParamsStruct.h` | `configParamStruct` — 32 run-time parameters |
| `detectorStructs.h` | `detectorSpecStruct` (19+2 fields), `vecField2D` (E-field), `scalarField3D` (phi) |
| `trailLogStruct.h` | `trailLog` — per-charge-element trajectory (10 vectors: x/y/z + q for e/h) |
| `phiLogStruct.h` | `phiLog` — weighting potential along trails (EPhi, HPhi) |

### Build Command

```bash
bash compile_all.sh
```

Or manually:
```bash
g++ -c src/Detector.cpp src/CTSI.cpp src/main.cpp src/Event.cpp src/MyRandom.cpp src/Preamplification.cpp src/Simulation.cpp
g++ -o ctsi.exe main.o CTSI.o Detector.o Event.o MyRandom.o Preamplification.o Simulation.o
chmod a+x ctsi.exe
```

---

## 15. Shiloh Variant

*Source: `modified-version-Shiloh/csti_v_5/CTSI_input_txt_files/detectorSpecL4p5.txt`*

An alternate detector geometry exists in the `modified-version-Shiloh/` directory with different parameters:

| Parameter | Standard | Shiloh |
|-----------|----------|--------|
| `L` (thickness) | 0.5 cm | **0.1 cm** |
| `W` (width) | 4 cm | **0.5 cm** |
| `NUM_ANODES` | 39 | 39 |
| `NUM_CATHODES` | 8 | **7** |
| `CATHODE_PITCH` | 5 mm | **5.5714 mm** |
| `CATHODE_WIDTH` | 4900 um | **5400 um** |
| `BIAS` | 1200 V/cm | **1000 V/cm** |
| `TAU_E` | 1e-5 s | **2e-5 s** |
| `TAU_H` | 1e-6 s | **1.5e-6 s** |

This represents a much thinner, narrower detector with slightly different material properties. All other parameters (mobilities, Fano factor, temperature, etc.) remain the same.

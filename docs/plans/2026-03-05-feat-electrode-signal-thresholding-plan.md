---
title: "feat: Electrode signal thresholding to reduce computation"
type: feat
status: active
date: 2026-03-05
---

# feat: Electrode signal thresholding to reduce computation

## Overview

Currently, CTSI computes induced-current waveforms for **all 39 anodes and all 8 cathodes** regardless of interaction position. In practice, only 2-5 neighboring anodes and 1-2 cathodes produce meaningful signals. This wastes ~85% of signal computation time.

This plan adds two layers of filtering:
1. **C++ simulator**: Limit which electrodes compute full signal chains (big speed win)
2. **Python toolkit**: Auto-select significant electrodes for display (no manual `--anodes` flags)

## Problem Statement

**Performance:** `Simulation.cpp:556-559` forces all electrodes to compute signals:
```cpp
for (int i = 0; i < (int) anCollectorOneHot.size(); i++)
    anCollectorOneHot[i] = 1;
for (int i = 0; i < (int) caCollectorOneHot.size(); i++)
    caCollectorOneHot[i] = 1;
```

This means `getTrailPhi()` samples the weighting potential grid for **every electrode** at **every time step** for **every charge element** — the dominant cost in signal computation.

**Usability:** Plotters show all electrodes by default, requiring manual `--anodes 18 19 20` flags. Users want automatic selection of relevant electrodes.

## Proposed Solution

### Phase 1: C++ — Proximity-based electrode window

Replace the "force all" logic with a configurable proximity window around the primary collecting electrode.

**How it works:**
1. After charge transport, the code already knows which electrode(s) collected charge via `anCollectorOneHot[i]`
2. Instead of forcing all to 1, expand only to ±N neighbors of each collecting electrode
3. N is configurable via a new `ctsi.cfg` parameter

**New config parameter:**
```
33_Neighbor_electrode_window    3
```

Meaning: compute signals for collecting anode ± 3 neighbors (7 anodes total instead of 39). For cathodes, with only 8 strips, compute all cathodes (or ± 1 neighbor).

**Code change in `Simulation.cpp`:**

Replace lines 556-559 with:
```cpp
// Expand collector set to include neighboring electrodes
// within configurable window for neighbor signal analysis.
int window = configParams.NEIGHBOR_WINDOW;

// For anodes: find min/max collector, expand by ±window
if (window >= 0) {
    int anMin = NUM_ANODES, anMax = -1;
    for (int i = 0; i < (int) anCollectorOneHot.size(); i++) {
        if (anCollectorOneHot[i] == 1) {
            anMin = std::min(anMin, i);
            anMax = std::max(anMax, i);
        }
    }
    if (anMax >= 0) {  // at least one collector found
        int lo = std::max(0, anMin - window);
        int hi = std::min((int)anCollectorOneHot.size() - 1, anMax + window);
        for (int i = lo; i <= hi; i++)
            anCollectorOneHot[i] = 1;
    }
}

// For cathodes: always compute all (only 8 strips, negligible cost)
for (int i = 0; i < (int) caCollectorOneHot.size(); i++)
    caCollectorOneHot[i] = 1;
```

**Backward compatibility:**
- `NEIGHBOR_WINDOW = -1` → compute ALL electrodes (current behavior)
- `NEIGHBOR_WINDOW = 0` → only collecting electrodes
- `NEIGHBOR_WINDOW = 3` → collecting ± 3 neighbors (recommended default)

**Expected speedup:** Signal computation scales linearly with electrode count. Going from 39 → ~7 anodes = **~5.5x faster** signal computation phase.

### Phase 2: Python — Auto-threshold for display

Add automatic electrode selection to all plotters based on peak signal amplitude.

**Threshold logic** (new utility function):
```python
def auto_select_electrodes(
    result: SimulationResult,
    threshold: float = 0.01,  # 1% of max signal
) -> tuple[list[int], list[int]]:
    """Return (anode_ids, cathode_ids) with significant signals.

    Selects electrodes whose peak absolute signal exceeds
    `threshold` fraction of the strongest electrode's peak.
    """
```

**Where it applies:**
- `plot_waveforms()`: replace default "all collectors" with auto-selected
- `plot_event()`: filter panels D & E to significant electrodes only
- `plot_zscan()`: auto-select electrodes if `--electrodes` not provided

**CLI changes:**
- Add `--auto-threshold` flag (default: ON) to plot commands
- Add `--threshold` float (default: 0.01) for the amplitude fraction
- Existing `--anodes`/`--cathodes`/`--electrodes` flags override auto-selection

## Technical Considerations

### Weighting potential proxy (alternative considered)

We considered evaluating weighting potential at the interaction point as a proxy for signal amplitude before computing full signal chains. **Rejected because:**
- The induced current depends on the *gradient* of the weighting potential along the charge path, not the value at one point
- A single-point check is unreliable — edge electrodes can have low weighting potential at the start but significant gradient
- Geometric proximity is simpler, faster (no extra lookups), and more predictable

### Edge cases

- **Multi-hit events:** Multiple interactions at different x-positions → union of all neighbor windows
- **Edge anodes (0 or 38):** Window clamps to valid range with `std::min`/`std::max`
- **Window = -1 sentinel:** Must be explicitly checked to preserve current behavior for validation runs

### Output format

No changes to the output file format. The parser already handles variable-length collector ID arrays. With fewer electrodes:
- `an_collector_id` shrinks (e.g., 7 instead of 39)
- `an_vs_time` has fewer waveforms
- `anode_id` (triggered) is unchanged — still only triggered electrodes

## Acceptance Criteria

- [x] New `NEIGHBOR_WINDOW` parameter in `configParamsStruct.h` and parsed from `ctsi.cfg`
- [x] `Simulation.cpp` uses proximity window instead of forcing all electrodes
- [x] `NEIGHBOR_WINDOW = -1` reproduces exact current behavior (all electrodes)
- [x] Default `ctsi.cfg` ships with `NEIGHBOR_WINDOW = 3`
- [x] `auto_select_electrodes()` utility in Python toolkit
- [x] All three plotters use auto-selection by default
- [x] `--anodes`/`--cathodes`/`--electrodes` CLI flags override auto-selection
- [ ] CI workflow runs successfully with new defaults
- [ ] Z-scan at 5 depths completes faster than before

## Files to Modify

### C++ changes
- `ctsi/src/configParamsStruct.h` — add `NEIGHBOR_WINDOW` field
- `ctsi/src/CTSI.cpp` — parse parameter 33 from config
- `ctsi/src/Simulation.cpp:556-559` — replace force-all with proximity window
- `ctsi/config/ctsi.cfg` — add line 33

### Python changes
- `ctsi/tools/ctsi_toolkit/plotters/__init__.py` — add `auto_select_electrodes()`
- `ctsi/tools/ctsi_toolkit/plotters/waveforms.py` — use auto-select as default
- `ctsi/tools/ctsi_toolkit/plotters/event_viz.py` — filter to significant electrodes
- `ctsi/tools/ctsi_toolkit/plotters/zscan.py` — auto-select if no `--electrodes`
- `ctsi/tools/ctsi_toolkit/cli.py` — add `--auto-threshold` and `--threshold` flags

## References

- Forcing logic: `ctsi/src/Simulation.cpp:556-559`
- Trigger thresholds: `ctsi/src/Simulation.cpp:1207-1238`
- Collector ID construction: `ctsi/src/Simulation.cpp:1188-1198`
- Waveform plotter selection: `ctsi/tools/ctsi_toolkit/plotters/waveforms.py:50-70`
- Z-scan electrode parsing: `ctsi/tools/ctsi_toolkit/plotters/zscan.py:49-57`
- Config parsing: `ctsi/src/CTSI.cpp:151-340`

"""Parse interactiveOut.txt into a typed SimulationResult dataclass.

Layout follows the authoritative MATLAB parser ``interactiveOut.m`` exactly.
Ported from the Python version in ``plot_simulation.py:13-244``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class TrailLog:
    """Position and charge history for one charge element."""
    xPosE: np.ndarray
    yPosE: np.ndarray
    zPosE: np.ndarray
    qETrapped: np.ndarray
    qEMobile: np.ndarray
    xPosH: np.ndarray
    yPosH: np.ndarray
    zPosH: np.ndarray
    qHTrapped: np.ndarray
    qHMobile: np.ndarray


@dataclass
class SimulationResult:
    """All data produced by a single ctsi.exe interactive-mode run."""

    # Time axes
    time_vec: np.ndarray                       # (time_len,)
    time_vec_preamp: np.ndarray                # (preamp_len,)

    # Collector IDs (electrodes that collected charge)
    an_collector_id: np.ndarray                # int array
    ca_collector_id: np.ndarray                # int array

    # Triggered electrode info
    anode_id: np.ndarray                       # int array
    cathode_id: np.ndarray                     # int array
    an_energy: np.ndarray
    ca_energy: np.ndarray
    an_trigger_time: np.ndarray
    ca_trigger_time: np.ndarray
    noisy_an_energy: np.ndarray
    noisy_ca_energy: np.ndarray
    noisy_an_trigger_time: np.ndarray
    noisy_ca_trigger_time: np.ndarray

    # Induced-current waveforms (one per collector electrode)
    an_vs_time: List[np.ndarray] = field(default_factory=list)
    ca_vs_time: List[np.ndarray] = field(default_factory=list)

    # Event trails: trails[interaction][charge_element] → TrailLog
    event_trails: List[List[TrailLog]] = field(default_factory=list)

    # Trail sizes
    trail_size_e: List[List[int]] = field(default_factory=list)
    trail_size_h: List[List[int]] = field(default_factory=list)

    # Preamp-domain signals (one per collector electrode)
    an_vs_time_reg: List[np.ndarray] = field(default_factory=list)
    ca_vs_time_reg: List[np.ndarray] = field(default_factory=list)
    noisy_an_vs_time_reg: List[np.ndarray] = field(default_factory=list)
    noisy_ca_vs_time_reg: List[np.ndarray] = field(default_factory=list)
    an_vs_time_preamp: List[np.ndarray] = field(default_factory=list)
    ca_vs_time_preamp: List[np.ndarray] = field(default_factory=list)
    noisy_an_vs_time_preamp: List[np.ndarray] = field(default_factory=list)
    noisy_ca_vs_time_preamp: List[np.ndarray] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

_TRAIL_FIELDS = [
    "xPosE", "yPosE", "zPosE", "qETrapped", "qEMobile",
    "xPosH", "yPosH", "zPosH", "qHTrapped", "qHMobile",
]


def parse_interactive_output(
    filename: str | Path = "interactiveOut.txt",
    *,
    mode: str = "full",
) -> SimulationResult:
    """Parse interactiveOut.txt and return a :class:`SimulationResult`.

    Parameters
    ----------
    filename : str or Path
        Path to the output file produced by ``ctsi.exe i``.
    mode : str
        ``"full"``  — parse everything (trails, phi, preamp signals).
        ``"waveforms_only"`` — parse up through induced-current waveforms,
        skip trails / phi / preamp for speed.
    """
    data = np.loadtxt(filename)
    idx = 0

    # ------------------------------------------------------------------
    # Time vector
    # ------------------------------------------------------------------
    time_len = int(data[idx]); idx += 1
    time_vec = data[idx:idx + time_len]; idx += time_len

    # ------------------------------------------------------------------
    # Collector IDs
    # ------------------------------------------------------------------
    num_an_collected = int(data[idx]); idx += 1
    an_collector_id = data[idx:idx + num_an_collected].astype(int) if num_an_collected > 0 else np.array([], dtype=int)
    idx += num_an_collected

    num_ca_collected = int(data[idx]); idx += 1
    ca_collector_id = data[idx:idx + num_ca_collected].astype(int) if num_ca_collected > 0 else np.array([], dtype=int)
    idx += num_ca_collected

    # ------------------------------------------------------------------
    # Anode trigger section
    # ------------------------------------------------------------------
    num_an_trig = int(data[idx]); idx += 1

    def _read_vec(n):
        nonlocal idx
        v = data[idx:idx + n] if n > 0 else np.array([])
        idx += n
        return v

    def _read_int_vec(n):
        nonlocal idx
        v = data[idx:idx + n].astype(int) if n > 0 else np.array([], dtype=int)
        idx += n
        return v

    anode_id = _read_int_vec(num_an_trig)
    an_energy = _read_vec(num_an_trig)
    an_trigger_time = _read_vec(num_an_trig)
    noisy_an_energy = _read_vec(num_an_trig)
    noisy_an_trigger_time = _read_vec(num_an_trig)

    # anVsTime — one waveform per collected anode
    num_an_signals = int(data[idx]); idx += 1
    an_vs_time = []
    for _ in range(num_an_signals):
        an_vs_time.append(data[idx:idx + time_len]); idx += time_len

    # ------------------------------------------------------------------
    # Cathode trigger section
    # ------------------------------------------------------------------
    num_ca_trig = int(data[idx]); idx += 1
    cathode_id = _read_int_vec(num_ca_trig)
    ca_energy = _read_vec(num_ca_trig)
    ca_trigger_time = _read_vec(num_ca_trig)
    noisy_ca_energy = _read_vec(num_ca_trig)
    noisy_ca_trigger_time = _read_vec(num_ca_trig)

    # caVsTime — one waveform per collected cathode
    num_ca_signals = int(data[idx]); idx += 1
    ca_vs_time = []
    for _ in range(num_ca_signals):
        ca_vs_time.append(data[idx:idx + time_len]); idx += time_len

    # Early exit for waveforms-only mode
    if mode == "waveforms_only":
        return SimulationResult(
            time_vec=time_vec,
            time_vec_preamp=np.array([]),
            an_collector_id=an_collector_id,
            ca_collector_id=ca_collector_id,
            anode_id=anode_id,
            cathode_id=cathode_id,
            an_energy=an_energy,
            ca_energy=ca_energy,
            an_trigger_time=an_trigger_time,
            ca_trigger_time=ca_trigger_time,
            noisy_an_energy=noisy_an_energy,
            noisy_ca_energy=noisy_ca_energy,
            noisy_an_trigger_time=noisy_an_trigger_time,
            noisy_ca_trigger_time=noisy_ca_trigger_time,
            an_vs_time=an_vs_time,
            ca_vs_time=ca_vs_time,
        )

    # ------------------------------------------------------------------
    # Event trails
    # ------------------------------------------------------------------
    num_fields_in_trail = int(data[idx]); idx += 1
    num_intrxn = int(data[idx]); idx += 1

    event_trails: List[List[TrailLog]] = []
    for _ in range(num_intrxn):
        num_charge_elems = int(data[idx]); idx += 1
        intrxn_trails: List[TrailLog] = []
        for _ in range(num_charge_elems):
            fields = {}
            for fname in _TRAIL_FIELDS:
                length = int(data[idx]); idx += 1
                fields[fname] = data[idx:idx + length] if length > 0 else np.array([])
                idx += length
            intrxn_trails.append(TrailLog(**fields))
        event_trails.append(intrxn_trails)

    # ------------------------------------------------------------------
    # Phi sections (anode then cathode) — skip contents, just advance idx
    # ------------------------------------------------------------------
    # Anode phi
    num_anodes_phi = int(data[idx]); idx += 1
    for _ in range(num_anodes_phi):
        num_intrxn_phi = int(data[idx]); idx += 1
        for _ in range(num_intrxn_phi):
            num_ce = int(data[idx]); idx += 1
            for _ in range(num_ce):
                length = int(data[idx]); idx += 1 + length  # EPhi
                length = int(data[idx]); idx += 1 + length  # HPhi

    # Cathode phi
    num_cathodes_phi = int(data[idx]); idx += 1
    for _ in range(num_cathodes_phi):
        num_intrxn_phi = int(data[idx]); idx += 1
        for _ in range(num_intrxn_phi):
            num_ce = int(data[idx]); idx += 1
            for _ in range(num_ce):
                length = int(data[idx]); idx += 1 + length  # EPhi
                length = int(data[idx]); idx += 1 + length  # HPhi

    # ------------------------------------------------------------------
    # Trail sizes (electron then hole)
    # ------------------------------------------------------------------
    trail_size_e: List[List[int]] = []
    num_intrxn_sizes = int(data[idx]); idx += 1
    for _ in range(num_intrxn_sizes):
        num_ce = int(data[idx]); idx += 1
        sizes = []
        for _ in range(num_ce):
            sizes.append(int(data[idx])); idx += 1
        trail_size_e.append(sizes)

    trail_size_h: List[List[int]] = []
    num_intrxn_sizes = int(data[idx]); idx += 1
    for _ in range(num_intrxn_sizes):
        num_ce = int(data[idx]); idx += 1
        sizes = []
        for _ in range(num_ce):
            sizes.append(int(data[idx])); idx += 1
        trail_size_h.append(sizes)

    # ------------------------------------------------------------------
    # Preamp time vector
    # ------------------------------------------------------------------
    preamp_len = int(data[idx]); idx += 1
    time_vec_preamp = data[idx:idx + preamp_len]; idx += preamp_len

    # ------------------------------------------------------------------
    # Regularised and preamp signals
    # ------------------------------------------------------------------
    def _read_signals(count):
        nonlocal idx
        out = []
        for _ in range(count):
            out.append(data[idx:idx + preamp_len]); idx += preamp_len
        return out

    an_vs_time_reg = _read_signals(num_an_collected)
    ca_vs_time_reg = _read_signals(num_ca_collected)
    noisy_an_vs_time_reg = _read_signals(num_an_collected)
    noisy_ca_vs_time_reg = _read_signals(num_ca_collected)
    an_vs_time_preamp = _read_signals(num_an_collected)
    ca_vs_time_preamp = _read_signals(num_ca_collected)
    noisy_an_vs_time_preamp = _read_signals(num_an_collected)
    noisy_ca_vs_time_preamp = _read_signals(num_ca_collected)

    return SimulationResult(
        time_vec=time_vec,
        time_vec_preamp=time_vec_preamp,
        an_collector_id=an_collector_id,
        ca_collector_id=ca_collector_id,
        anode_id=anode_id,
        cathode_id=cathode_id,
        an_energy=an_energy,
        ca_energy=ca_energy,
        an_trigger_time=an_trigger_time,
        ca_trigger_time=ca_trigger_time,
        noisy_an_energy=noisy_an_energy,
        noisy_ca_energy=noisy_ca_energy,
        noisy_an_trigger_time=noisy_an_trigger_time,
        noisy_ca_trigger_time=noisy_ca_trigger_time,
        an_vs_time=an_vs_time,
        ca_vs_time=ca_vs_time,
        event_trails=event_trails,
        trail_size_e=trail_size_e,
        trail_size_h=trail_size_h,
        an_vs_time_reg=an_vs_time_reg,
        ca_vs_time_reg=ca_vs_time_reg,
        noisy_an_vs_time_reg=noisy_an_vs_time_reg,
        noisy_ca_vs_time_reg=noisy_ca_vs_time_reg,
        an_vs_time_preamp=an_vs_time_preamp,
        ca_vs_time_preamp=ca_vs_time_preamp,
        noisy_an_vs_time_preamp=noisy_an_vs_time_preamp,
        noisy_ca_vs_time_preamp=noisy_ca_vs_time_preamp,
    )

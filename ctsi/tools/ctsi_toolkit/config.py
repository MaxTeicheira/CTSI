"""Parse detectorSpec.txt and ctsi.cfg into typed dataclasses."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# DetectorSpec — parsed from detectorSpec.txt
# ---------------------------------------------------------------------------

@dataclass
class DetectorSpec:
    """Detector geometry and material properties from detectorSpec.txt."""

    L: float               # Thickness (cm)
    W: float               # Width (cm)
    NUM_ANODES: int
    NUM_CATHODES: int
    ANODE_PITCH: float     # mm
    ANODE_WIDTH: float     # um
    CATHODE_PITCH: float   # mm
    CATHODE_WIDTH: float   # um
    BIAS: float            # V/cm
    CZT_W_FACTOR: float   # eV
    DIFFUSION_E: float     # cm^2/s
    DIFFUSION_H: float     # cm^2/s
    MU_E: float            # cm^2/Vs
    MU_H: float            # cm^2/Vs
    TAU_E: float           # s
    TAU_H: float           # s
    FANO: float
    EPSILON_R: float
    TEMP: float            # K

    # Derived convenience properties
    @property
    def L_mm(self) -> float:
        return self.L * 10.0

    @property
    def W_mm(self) -> float:
        return self.W * 10.0


# Type map for parsing — int fields vs float fields
_INT_FIELDS = {"NUM_ANODES", "NUM_CATHODES"}

_FIELD_TYPES = {
    "L": float, "W": float, "NUM_ANODES": int, "NUM_CATHODES": int,
    "ANODE_PITCH": float, "ANODE_WIDTH": float,
    "CATHODE_PITCH": float, "CATHODE_WIDTH": float,
    "BIAS": float, "CZT_W_FACTOR": float,
    "DIFFUSION_E": float, "DIFFUSION_H": float,
    "MU_E": float, "MU_H": float,
    "TAU_E": float, "TAU_H": float,
    "FANO": float, "EPSILON_R": float, "TEMP": float,
}


def load_detector_spec(path: str | Path = "config/detectorSpec.txt") -> DetectorSpec:
    """Parse ``KEY value // comment`` format into a DetectorSpec."""
    vals: dict = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("//") or line == "<":
                continue
            # Strip trailing comment
            if "//" in line:
                line = line[:line.index("//")]
            parts = line.split()
            if len(parts) < 2:
                continue
            key, raw_val = parts[0], parts[1]
            if key in _FIELD_TYPES:
                vals[key] = _FIELD_TYPES[key](float(raw_val))
    return DetectorSpec(**vals)


# ---------------------------------------------------------------------------
# CtsiConfig — parsed from ctsi.cfg
# ---------------------------------------------------------------------------

@dataclass
class CtsiConfig:
    """Runtime configuration from ctsi.cfg."""

    efield_path: str = ""
    weighting_potential_source: str = "n"
    anode_phi_path: str = ""
    cathode_phi_path: str = ""
    gray_output: str = ""
    detector_spec_path: str = ""
    preamp_spec_path: str = ""
    list_output: str = "listOut.txt"
    random_seed: int = 65539
    num_charge_elements: int = 50
    max_sim_steps: int = 10000
    small_pixel_distance: float = 0.445
    anode_trigger_threshold: float = 0.05
    cathode_trigger_threshold: float = -0.04
    fwhm_preamp_noise: float = 0.0
    efield_x_grid: float = 0.01
    efield_y_grid: float = 1.0
    efield_z_grid: float = 0.01
    anode_phi_x_grid: float = 0.01
    anode_phi_y_grid: float = 0.01
    anode_phi_z_grid: float = 0.01
    cathode_phi_x_grid: float = 0.01
    cathode_phi_y_grid: float = 0.01
    cathode_phi_z_grid: float = 0.01
    event_x_scale: float = 1.0
    event_y_scale: float = 1.0
    event_z_scale: float = -1.0
    event_x_offset: float = 0.0
    event_y_offset: float = 0.0
    event_z_offset: float = 0.5
    output_mode: str = "r"
    caution: str = "false"
    neighbor_electrode_window: int = 3


# Map ctsi.cfg line numbers to CtsiConfig field names
_CFG_MAP = {
    1:  "efield_path",
    2:  "weighting_potential_source",
    3:  "anode_phi_path",
    4:  "cathode_phi_path",
    5:  "gray_output",
    6:  "detector_spec_path",
    7:  "preamp_spec_path",
    8:  "list_output",
    9:  "random_seed",
    10: "num_charge_elements",
    11: "max_sim_steps",
    12: "small_pixel_distance",
    13: "anode_trigger_threshold",
    14: "cathode_trigger_threshold",
    15: "fwhm_preamp_noise",
    16: "efield_x_grid",
    17: "efield_y_grid",
    18: "efield_z_grid",
    19: "anode_phi_x_grid",
    20: "anode_phi_y_grid",
    21: "anode_phi_z_grid",
    22: "cathode_phi_x_grid",
    23: "cathode_phi_y_grid",
    24: "cathode_phi_z_grid",
    25: "event_x_scale",
    26: "event_y_scale",
    27: "event_z_scale",
    28: "event_x_offset",
    29: "event_y_offset",
    30: "event_z_offset",
    31: "output_mode",
    32: "caution",
    33: "neighbor_electrode_window",
}

_CFG_FLOAT_FIELDS = {
    "small_pixel_distance", "anode_trigger_threshold", "cathode_trigger_threshold",
    "fwhm_preamp_noise",
    "efield_x_grid", "efield_y_grid", "efield_z_grid",
    "anode_phi_x_grid", "anode_phi_y_grid", "anode_phi_z_grid",
    "cathode_phi_x_grid", "cathode_phi_y_grid", "cathode_phi_z_grid",
    "event_x_scale", "event_y_scale", "event_z_scale",
    "event_x_offset", "event_y_offset", "event_z_offset",
}

_CFG_INT_FIELDS = {"random_seed", "num_charge_elements", "max_sim_steps", "neighbor_electrode_window"}


def load_ctsi_config(path: str | Path = "config/ctsi.cfg") -> CtsiConfig:
    """Parse ``N_Label<tab>value`` format into a CtsiConfig."""
    cfg = CtsiConfig()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)  # split on any whitespace, max 2 parts
            if len(parts) < 2:
                continue
            label = parts[0]
            value = parts[1].strip()
            # Extract line number from label (e.g. "1_E_field" → 1)
            try:
                line_num = int(label.split("_", 1)[0])
            except ValueError:
                continue
            field_name = _CFG_MAP.get(line_num)
            if field_name is None:
                continue
            if field_name in _CFG_FLOAT_FIELDS:
                setattr(cfg, field_name, float(value))
            elif field_name in _CFG_INT_FIELDS:
                setattr(cfg, field_name, int(value))
            else:
                setattr(cfg, field_name, value)
    return cfg


def write_ctsi_config(cfg: CtsiConfig, path: str | Path = "config/ctsi.cfg") -> None:
    """Write a CtsiConfig back to disk in the expected format."""
    # Reverse map: field_name → (line_num, label_suffix)
    _LABELS = {
        1: "E_field", 2: "Weighting_potential_source",
        3: "Anode_weighting_potential", 4: "Cathode_weighting_potential",
        5: "GRAY_output", 6: "Detector_specification",
        7: "Preamplifier_specification", 8: "List-mode_output_name",
        9: "Random_seed", 10: "Number_of_charge_elements",
        11: "Max_num_sim_steps", 12: "Small_pixel_distance",
        13: "Anode_trigger_threshold", 14: "Cathode_trigger_threshold",
        15: "FWHM_preamplifier_noise", 16: "E_field_x_grid_space",
        17: "E_field_y_grid_space", 18: "E_field_z_grid_space",
        19: "Anode_phi_x_grid_space", 20: "Anode_phi_y_grid_space",
        21: "Anode_phi_z_grid_space", 22: "Cathode_phi_x_grid_space",
        23: "Cathode_phi_y_grid_space", 24: "Cathode_phi_z_grid_space",
        25: "Event_x_pos_scale_factor", 26: "Event_y_pos_scale_factor",
        27: "Event_z_pos_scale_factor", 28: "Event_x_pos_offset",
        29: "Event_y_pos_offset", 30: "Event_z_pos_offset",
        31: "Output_mode", 32: "Caution",
        33: "Neighbor_electrode_window",
    }
    with open(path, "w") as f:
        for line_num in sorted(_CFG_MAP.keys()):
            field_name = _CFG_MAP[line_num]
            value = getattr(cfg, field_name)
            label = f"{line_num}_{_LABELS[line_num]}"
            f.write(f"{label}\t\t{value}\n")


# ---------------------------------------------------------------------------
# Event file creation
# ---------------------------------------------------------------------------

def create_event_file(
    x: float,
    y: float,
    z: float,
    energy: float,
    cfg: CtsiConfig,
    path: str | Path = "single_event.txt",
    n_interactions: int = 1,
) -> Path:
    """Write a single-event file with proper z-coordinate transform.

    The z value is in *detector* coordinates (0 = cathode, L = anode).
    The transform z_event = z_scale * z_detector + z_offset is applied
    so that ctsi.exe reads the correct position.

    Parameters
    ----------
    x, y, z : float
        Interaction position in detector coordinates (cm).
    energy : float
        Gamma-ray energy in MeV.
    cfg : CtsiConfig
        Used to compute the z transform (scale + offset).
    path : str or Path
        Output file path.
    n_interactions : int
        Number of interactions (normally 1 for single-gamma).

    Returns
    -------
    Path to the written event file.
    """
    # Inverse transform: z_event such that z_scale * z_event + z_offset = z_detector
    # => z_event = (z_detector - z_offset) / z_scale
    z_event = (z - cfg.event_z_offset) / cfg.event_z_scale

    path = Path(path)
    with open(path, "w") as f:
        # Format: type interactions ? time energy x y z ? ?
        f.write(f"3\t{n_interactions}\t0\t1.0e-9\t{energy}\t{x}\t{y}\t{z_event}\t1\t1\n")
    return path

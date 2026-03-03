"""Run ctsi.exe simulations and return parsed results.

Provides two modes:
  - ``run_single_event``: one interaction, pipes stdin to ``ctsi.exe i``
  - ``run_zscan``: loops over z positions, returns dict of z → SimulationResult
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .config import CtsiConfig, DetectorSpec, load_ctsi_config, load_detector_spec
from .parser import SimulationResult, parse_interactive_output


def run_single_event(
    x: float,
    y: float,
    z: float,
    energy: float,
    *,
    exe: str = "./ctsi.exe",
    output_file: str = "output/interactiveOut.txt",
    timeout: int = 300,
    mode: str = "full",
) -> SimulationResult:
    """Run a single event via interactive mode and return the parsed result.

    Parameters
    ----------
    x, y, z : float
        Interaction position in *detector* coordinates (cm).
        z=0 is cathode, z=L is anode.
    energy : float
        Gamma-ray energy in MeV.
    exe : str
        Path to the ctsi executable.
    output_file : str
        Where ctsi.exe writes its interactive output.
    timeout : int
        Max seconds to wait for the simulation.
    mode : str
        Parser mode — ``"full"`` or ``"waveforms_only"``.

    Returns
    -------
    SimulationResult
    """
    stdin_text = f"{x}\n{y}\n{z}\n{energy}\nn\n"

    proc = subprocess.run(
        [exe, "i"],
        input=stdin_text,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"ctsi.exe exited with code {proc.returncode}\n"
            f"stderr: {proc.stderr[:500]}"
        )

    return parse_interactive_output(output_file, mode=mode)


def run_zscan(
    x: float,
    y: float,
    z_positions: List[float] | np.ndarray,
    energy: float,
    *,
    exe: str = "./ctsi.exe",
    output_file: str = "output/interactiveOut.txt",
    save_dir: Optional[str] = None,
    timeout: int = 300,
    mode: str = "waveforms_only",
    verbose: bool = True,
) -> Dict[float, SimulationResult]:
    """Run events at multiple depths and return a z → SimulationResult dict.

    Parameters
    ----------
    x, y : float
        Fixed lateral position (cm).
    z_positions : array-like
        Depths to simulate in detector coordinates (cm).
    energy : float
        Gamma-ray energy in MeV.
    exe : str
        Path to ctsi executable.
    output_file : str
        Intermediate output file that ctsi.exe writes each run.
    save_dir : str, optional
        If given, copy each interactiveOut.txt into ``save_dir/z_<val>.txt``
        for reproducibility.
    timeout : int
        Per-event timeout in seconds.
    mode : str
        Parser mode (``"waveforms_only"`` is faster for z-scan plots).
    verbose : bool
        Print progress to stdout.

    Returns
    -------
    dict mapping z_position (float) → SimulationResult
    """
    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)

    results: Dict[float, SimulationResult] = {}

    for i, z in enumerate(z_positions):
        if verbose:
            print(f"[{i + 1}/{len(z_positions)}] z = {z:.3f} cm ({z * 10:.1f} mm)...")

        result = run_single_event(
            x, y, z, energy,
            exe=exe, output_file=output_file, timeout=timeout, mode=mode,
        )
        results[float(z)] = result

        if save_dir:
            dest = save_path / f"z_{z:.4f}.txt"
            shutil.copy2(output_file, dest)

        if verbose:
            # Report peak anode signal
            if result.an_vs_time:
                peak = max(np.max(np.abs(wf)) for wf in result.an_vs_time)
                print(f"   Peak anode signal: {peak:.4g}")

    if verbose:
        print("Done.")

    return results

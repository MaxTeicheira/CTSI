"""Depth-dependent 4-panel waveform plot with jet colorbar.

Replaces: the plotting section of run_and_plot_zscan.py
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable, jet
from matplotlib.colors import Normalize

from ..parser import SimulationResult


def plot_zscan(
    z_results: Dict[float, SimulationResult],
    electrodes: Optional[List[str]] = None,
    *,
    output: Optional[str] = None,
    title: Optional[str] = None,
    time_unit: str = "ns",
    show: bool = False,
) -> plt.Figure:
    """Plot a multi-depth waveform comparison with jet colormap.

    Parameters
    ----------
    z_results : dict
        Mapping of z_position (cm) → SimulationResult.
    electrodes : list of str, optional
        Electrode labels like ``["A18", "A19", "A20", "C1"]``.
        Format: ``A<n>`` for anode, ``C<n>`` for cathode.
        Defaults to ``["A18", "A19", "A20", "C1"]``.
    output : str, optional
        Save figure to this path.
    title : str, optional
        Custom suptitle.
    time_unit : str
        ``"ns"`` (default) or ``"us"``.
    show : bool
        Call ``plt.show()``.
    """
    if electrodes is None:
        electrodes = ["A18", "A19", "A20", "C1"]

    # Parse electrode specs
    parsed: List[Tuple[str, int]] = []  # ("anode"|"cathode", id)
    for e in electrodes:
        if e.upper().startswith("A"):
            parsed.append(("anode", int(e[1:])))
        elif e.upper().startswith("C"):
            parsed.append(("cathode", int(e[1:])))
        else:
            raise ValueError(f"Electrode format must be A<n> or C<n>, got: {e}")

    z_positions = sorted(z_results.keys())
    z_mm = [z * 10 for z in z_positions]

    # Time conversion
    scale = {"ns": 1e9, "us": 1e6}[time_unit]
    unit_label = {"ns": "ns", "us": "\u00b5s"}[time_unit]

    # Color normalization
    norm = Normalize(vmin=min(z_mm), vmax=max(z_mm))
    sm = ScalarMappable(cmap=jet, norm=norm)

    # Layout: ceil(n/2) rows x 2 cols
    n = len(parsed)
    nrows = (n + 1) // 2
    ncols = min(n, 2)

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 5 * nrows), sharex=True)
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    axes_flat = axes.flat

    for ax, (etype, eid) in zip(axes_flat, parsed):
        for z in z_positions:
            res = z_results[z]
            time = res.time_vec * scale
            color = jet(norm(z * 10))

            if etype == "anode":
                ids = list(res.an_collector_id)
                if eid in ids:
                    wf = res.an_vs_time[ids.index(eid)]
                    ax.plot(time, wf, color=color, label=f"z = {z*10:.1f} mm")
                label = f"Anode {eid}"
            else:
                ids = list(res.ca_collector_id)
                if eid in ids:
                    wf = res.ca_vs_time[ids.index(eid)]
                    ax.plot(time, wf, color=color, label=f"z = {z*10:.1f} mm")
                label = f"Cathode {eid}"

            ax.set_title(label, fontsize=13, fontweight="bold")
            ax.set_ylabel("Induced Current", fontsize=10)
            ax.grid(True, alpha=0.3)

    # X-axis label on bottom row only
    for ax in axes[-1]:
        ax.set_xlabel(f"Time ({unit_label})", fontsize=11)

    # Hide unused axes
    for i in range(n, nrows * ncols):
        axes_flat[i].set_visible(False)

    fig.suptitle(
        title or "Z-Scan: Waveforms at Multiple Depths",
        fontsize=14, fontweight="bold",
    )

    plt.tight_layout(rect=[0, 0, 0.88, 0.93])

    # Colorbar on right
    cbar_ax = fig.add_axes([0.90, 0.12, 0.02, 0.75])
    fig.colorbar(sm, cax=cbar_ax, label="Depth z (mm)")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved to {output}")
    if show:
        plt.show()

    return fig

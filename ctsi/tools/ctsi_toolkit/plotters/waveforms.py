"""Plot selected electrode waveforms from a single event.

Replaces: plot_single_event.py
"""

from __future__ import annotations

from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

from ..parser import SimulationResult


def plot_waveforms(
    result: SimulationResult,
    anodes: Optional[List[int]] = None,
    cathodes: Optional[List[int]] = None,
    *,
    output: Optional[str] = None,
    title: Optional[str] = None,
    time_unit: str = "ns",
    show: bool = False,
) -> plt.Figure:
    """Plot induced-current waveforms for selected electrodes.

    Parameters
    ----------
    result : SimulationResult
        Parsed simulation output.
    anodes : list of int, optional
        Anode IDs to plot (e.g. [18, 19, 20]). If None, plot all collectors.
    cathodes : list of int, optional
        Cathode IDs to plot (e.g. [1]). If None, plot all collectors.
    output : str, optional
        Save figure to this path.  If None, figure is returned without saving.
    title : str, optional
        Custom figure title.
    time_unit : str
        ``"ns"`` (default) or ``"us"`` for microseconds.
    show : bool
        Call ``plt.show()`` after plotting.
    """
    # Time conversion
    scale = {"ns": 1e9, "us": 1e6}[time_unit]
    unit_label = {"ns": "ns", "us": "\u00b5s"}[time_unit]
    time = result.time_vec * scale

    # Default to all collectors if no selection given
    an_ids = list(result.an_collector_id)
    ca_ids = list(result.ca_collector_id)

    if anodes is None:
        anodes = an_ids
    if cathodes is None:
        cathodes = ca_ids

    fig, ax = plt.subplots(figsize=(12, 7))

    for aid in anodes:
        if aid in an_ids:
            i = an_ids.index(aid)
            ax.plot(time, result.an_vs_time[i], linewidth=1.5, label=f"Anode {aid}")

    for cid in cathodes:
        if cid in ca_ids:
            i = ca_ids.index(cid)
            ax.plot(time, result.ca_vs_time[i], linewidth=1.5, linestyle="--",
                    label=f"Cathode {cid}")

    ax.set_xlabel(f"Time ({unit_label})", fontsize=12)
    ax.set_ylabel("Induced Current", fontsize=12)
    ax.set_title(title or "Single Event Waveforms", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved to {output}")
    if show:
        plt.show()

    return fig

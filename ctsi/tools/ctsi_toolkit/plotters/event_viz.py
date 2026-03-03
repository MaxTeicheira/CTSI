"""Full event visualization: trajectories + charge + signals.

Replaces: plot_event.py, plot_simulation.py, visualize_results.py
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ..config import DetectorSpec
from ..parser import SimulationResult


# Consistent color scheme
_E_COLOR = "#2166ac"   # blue — electrons
_H_COLOR = "#b2182b"   # red — holes


def plot_event(
    result: SimulationResult,
    spec: DetectorSpec,
    *,
    output: Optional[str] = None,
    title: Optional[str] = None,
    show: bool = False,
) -> plt.Figure:
    """Create a 6-panel event visualisation.

    Panels:
        (A) X-Z trajectories
        (B) Z vs time
        (C) Mobile charge vs time
        (D) Anode induced signals
        (E) Cathode induced signals
        (F) Summary statistics

    Parameters
    ----------
    result : SimulationResult
        Parsed simulation output (must be full mode with trails).
    spec : DetectorSpec
        Detector geometry for axis limits and summary text.
    output : str, optional
        Save figure to this path.
    title : str, optional
        Custom suptitle.
    show : bool
        Call ``plt.show()``.
    """
    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    time_us = result.time_vec * 1e6
    L_mm = spec.L_mm

    # ------------------------------------------------------------------
    # (A) X-Z trajectories
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    if result.event_trails:
        for i, trail in enumerate(result.event_trails[0]):
            alpha = 0.7 if i < 5 else 0.3
            lw = 1.5 if i < 5 else 0.8
            if len(trail.zPosE) > 0:
                ax.plot(trail.xPosE * 10, trail.zPosE * 10,
                        color=_E_COLOR, alpha=alpha, linewidth=lw)
            if len(trail.zPosH) > 0:
                ax.plot(trail.xPosH * 10, trail.zPosH * 10,
                        color=_H_COLOR, alpha=alpha, linewidth=lw, linestyle="--")

    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
    ax.axhline(y=L_mm, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("X Position (mm)")
    ax.set_ylabel("Z Position (mm)")
    ax.set_title("(A) Charge Trajectories (X-Z)")
    ax.set_ylim(-0.3, L_mm + 0.3)
    ax.plot([], [], color=_E_COLOR, lw=2, label="Electrons")
    ax.plot([], [], color=_H_COLOR, lw=2, ls="--", label="Holes")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (B) Z vs time
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    if result.event_trails:
        for i, trail in enumerate(result.event_trails[0][:10]):
            alpha = 0.8 - min(i, 5) * 0.1
            if len(trail.zPosE) > 0:
                t = time_us[:len(trail.zPosE)]
                ax.plot(t, trail.zPosE * 10, color=_E_COLOR, alpha=alpha, lw=1.5)
            if len(trail.zPosH) > 0:
                t = time_us[:len(trail.zPosH)]
                ax.plot(t, trail.zPosH * 10, color=_H_COLOR, alpha=alpha, lw=1.5, ls="--")

    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
    ax.axhline(y=L_mm, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Time (\u00b5s)")
    ax.set_ylabel("Z Position (mm)")
    ax.set_title("(B) Depth vs Time")
    ax.set_ylim(-0.3, L_mm + 0.3)
    ax.annotate("Anode", xy=(time_us[-1] * 0.8, L_mm), fontsize=9, color="gray", ha="center", va="bottom")
    ax.annotate("Cathode", xy=(time_us[-1] * 0.8, 0), fontsize=9, color="gray", ha="center", va="top")
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (C) Mobile charge vs time
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    if result.event_trails:
        total_e = np.zeros(len(time_us))
        total_h = np.zeros(len(time_us))
        for trail in result.event_trails[0]:
            if len(trail.qEMobile) > 0:
                total_e[:len(trail.qEMobile)] += trail.qEMobile
            if len(trail.qHMobile) > 0:
                total_h[:len(trail.qHMobile)] += trail.qHMobile

        ax.fill_between(time_us, 0, total_e, color=_E_COLOR, alpha=0.3)
        ax.fill_between(time_us, 0, total_h, color=_H_COLOR, alpha=0.3)
        ax.plot(time_us, total_e, color=_E_COLOR, lw=2, label="Electrons")
        ax.plot(time_us, total_h, color=_H_COLOR, lw=2, label="Holes")

    ax.set_xlabel("Time (\u00b5s)")
    ax.set_ylabel("Mobile Charge (a.u.)")
    ax.set_title("(C) Mobile Charge vs Time")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (D) Anode signals
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    for i, signal in enumerate(result.an_vs_time):
        an_id = int(result.an_collector_id[i]) if i < len(result.an_collector_id) else i
        ax.plot(time_us, signal, lw=1.5, label=f"Pixel {an_id}")
    ax.set_xlabel("Time (\u00b5s)")
    ax.set_ylabel("Induced Current (a.u.)")
    ax.set_title("(D) Anode Induced Signals")
    if result.an_vs_time:
        ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (E) Cathode signals
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    for i, signal in enumerate(result.ca_vs_time):
        ca_id = int(result.ca_collector_id[i]) if i < len(result.ca_collector_id) else i
        ax.plot(time_us, signal, lw=1.5, label=f"Strip {ca_id}")
    ax.set_xlabel("Time (\u00b5s)")
    ax.set_ylabel("Induced Current (a.u.)")
    ax.set_title("(E) Cathode Induced Signals")
    if result.ca_vs_time:
        ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # (F) Summary
    # ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")

    z_start = ""
    if result.event_trails and result.event_trails[0]:
        t0 = result.event_trails[0][0]
        if len(t0.zPosE) > 0:
            z_start = f"{t0.zPosE[0] * 10:.2f} mm"

    summary = (
        f"Simulation Summary\n"
        f"{'=' * 28}\n\n"
        f"Detector:  L = {L_mm:.1f} mm, W = {spec.W_mm:.0f} mm\n"
        f"Bias:      {spec.BIAS:.0f} V/cm\n\n"
        f"Time steps: {len(result.time_vec)}\n"
        f"Sim time:   {result.time_vec[-1] * 1e6:.2f} \u00b5s\n\n"
        f"Triggered anodes:   {list(result.anode_id)}\n"
        f"Triggered cathodes: {list(result.cathode_id)}\n\n"
        f"Anode energies:  {[f'{e:.4f}' for e in result.an_energy]}\n"
        f"Cathode energies: {[f'{e:.4f}' for e in result.ca_energy]}\n"
    )
    if z_start:
        summary += f"\nInteraction depth: {z_start}"

    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
            va="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle(
        title or "CZT Detector Simulation Results",
        fontsize=14, fontweight="bold", y=0.98,
    )

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved to {output}")
    if show:
        plt.show()

    return fig

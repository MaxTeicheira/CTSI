"""Argparse CLI entry point for the CTSI toolkit.

Usage::

    python -m ctsi_toolkit plot-waveforms --input interactiveOut.txt --anodes 18 19 20 --cathodes 1
    python -m ctsi_toolkit plot-zscan --electrodes A18 A19 A20 C1 --z-min 0.05 --z-max 0.45 --z-steps 5
    python -m ctsi_toolkit plot-event --input interactiveOut.txt
    python -m ctsi_toolkit run-event --x 0.01 --y 0.25 --z 0.25 --energy 0.662
    python -m ctsi_toolkit run-zscan --x 0.01 --y 0.25 --z-min 0.05 --z-max 0.45 --z-steps 5
"""

from __future__ import annotations

import argparse
import sys

import numpy as np


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add options shared across all subcommands."""
    parser.add_argument("--detector-spec", default="config/detectorSpec.txt",
                        help="Path to detectorSpec.txt")
    parser.add_argument("--config", default="config/ctsi.cfg",
                        help="Path to ctsi.cfg")
    parser.add_argument("--exe", default="./ctsi.exe",
                        help="Path to the ctsi executable")


# ── plot-waveforms ─────────────────────────────────────────────────────

def cmd_plot_waveforms(args: argparse.Namespace) -> None:
    from .parser import parse_interactive_output
    from .plotters.waveforms import plot_waveforms

    result = parse_interactive_output(args.input, mode="waveforms_only")
    plot_waveforms(
        result,
        anodes=args.anodes,
        cathodes=args.cathodes,
        output=args.output,
        title=args.title,
        show=args.show,
    )


# ── plot-zscan ─────────────────────────────────────────────────────────

def cmd_plot_zscan(args: argparse.Namespace) -> None:
    from .config import load_detector_spec
    from .plotters.zscan import plot_zscan
    from .runners import run_zscan

    z_positions = np.linspace(args.z_min, args.z_max, args.z_steps)
    z_results = run_zscan(
        args.x, args.y, z_positions, args.energy,
        exe=args.exe,
        save_dir=args.save_dir,
    )

    plot_zscan(
        z_results,
        electrodes=args.electrodes,
        output=args.output,
        title=args.title,
        show=args.show,
    )


# ── plot-event ─────────────────────────────────────────────────────────

def cmd_plot_event(args: argparse.Namespace) -> None:
    from .config import load_detector_spec
    from .parser import parse_interactive_output
    from .plotters.event_viz import plot_event

    spec = load_detector_spec(args.detector_spec)
    result = parse_interactive_output(args.input, mode="full")
    plot_event(
        result, spec,
        output=args.output,
        title=args.title,
        show=args.show,
    )


# ── run-event ──────────────────────────────────────────────────────────

def cmd_run_event(args: argparse.Namespace) -> None:
    from .runners import run_single_event

    result = run_single_event(
        args.x, args.y, args.z, args.energy,
        exe=args.exe,
    )

    print(f"Time steps:         {len(result.time_vec)}")
    print(f"Triggered anodes:   {list(result.anode_id)}")
    print(f"Triggered cathodes: {list(result.cathode_id)}")
    print(f"Anode energies:     {list(result.an_energy)}")
    print(f"Cathode energies:   {list(result.ca_energy)}")


# ── run-zscan ──────────────────────────────────────────────────────────

def cmd_run_zscan(args: argparse.Namespace) -> None:
    from .runners import run_zscan

    z_positions = np.linspace(args.z_min, args.z_max, args.z_steps)
    z_results = run_zscan(
        args.x, args.y, z_positions, args.energy,
        exe=args.exe,
        save_dir=args.save_dir,
    )

    print(f"\nCompleted {len(z_results)} positions.")
    for z in sorted(z_results):
        r = z_results[z]
        peak_an = max((np.max(np.abs(wf)) for wf in r.an_vs_time), default=0)
        print(f"  z={z*10:.1f} mm  anodes={list(r.an_collector_id)}  peak={peak_an:.4g}")


# ── main ───────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="ctsi_toolkit",
        description="CTSI CZT Detector Simulation Toolkit",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # -- plot-waveforms --
    p = sub.add_parser("plot-waveforms", help="Plot selected electrode waveforms")
    p.add_argument("--input", "-i", default="interactiveOut.txt",
                   help="Path to interactiveOut.txt")
    p.add_argument("--anodes", type=int, nargs="*", default=None,
                   help="Anode IDs to plot (e.g. 18 19 20)")
    p.add_argument("--cathodes", type=int, nargs="*", default=None,
                   help="Cathode IDs to plot (e.g. 1)")
    p.add_argument("--output", "-o", default="waveforms.png")
    p.add_argument("--title", default=None)
    p.add_argument("--show", action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_plot_waveforms)

    # -- plot-zscan --
    p = sub.add_parser("plot-zscan", help="Run z-scan and plot 4-panel figure")
    p.add_argument("--electrodes", nargs="*", default=["A18", "A19", "A20", "C1"],
                   help="Electrodes (e.g. A18 A19 C1)")
    p.add_argument("--x", type=float, default=0.01)
    p.add_argument("--y", type=float, default=0.25)
    p.add_argument("--energy", type=float, default=0.662)
    p.add_argument("--z-min", type=float, default=0.05)
    p.add_argument("--z-max", type=float, default=0.45)
    p.add_argument("--z-steps", type=int, default=5)
    p.add_argument("--save-dir", default=None,
                   help="Save intermediate outputs for reproducibility")
    p.add_argument("--output", "-o", default="zscan_4panel.png")
    p.add_argument("--title", default=None)
    p.add_argument("--show", action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_plot_zscan)

    # -- plot-event --
    p = sub.add_parser("plot-event", help="Full 6-panel event visualisation")
    p.add_argument("--input", "-i", default="interactiveOut.txt",
                   help="Path to interactiveOut.txt")
    p.add_argument("--output", "-o", default="event_viz.png")
    p.add_argument("--title", default=None)
    p.add_argument("--show", action="store_true")
    _add_common_args(p)
    p.set_defaults(func=cmd_plot_event)

    # -- run-event --
    p = sub.add_parser("run-event", help="Run a single simulation event")
    p.add_argument("--x", type=float, required=True)
    p.add_argument("--y", type=float, required=True)
    p.add_argument("--z", type=float, required=True)
    p.add_argument("--energy", type=float, required=True)
    _add_common_args(p)
    p.set_defaults(func=cmd_run_event)

    # -- run-zscan --
    p = sub.add_parser("run-zscan", help="Run z-scan sweep")
    p.add_argument("--x", type=float, default=0.01)
    p.add_argument("--y", type=float, default=0.25)
    p.add_argument("--energy", type=float, default=0.662)
    p.add_argument("--z-min", type=float, default=0.05)
    p.add_argument("--z-max", type=float, default=0.45)
    p.add_argument("--z-steps", type=int, default=5)
    p.add_argument("--save-dir", default=None)
    _add_common_args(p)
    p.set_defaults(func=cmd_run_zscan)

    args = parser.parse_args(argv)
    args.func(args)

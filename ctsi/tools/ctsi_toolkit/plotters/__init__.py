"""Plotting modules for CTSI simulation results."""

from __future__ import annotations

from typing import List, Tuple

import numpy as np


def auto_select_electrodes(
    result,
    threshold: float = 0.01,
) -> Tuple[List[int], List[int]]:
    """Return (anode_ids, cathode_ids) with significant signals.

    Selects electrodes whose peak absolute induced-current signal
    exceeds *threshold* fraction of the strongest electrode's peak.

    Parameters
    ----------
    result : SimulationResult
        Parsed simulation output.
    threshold : float
        Minimum fraction of the max peak signal to include (default 1%).

    Returns
    -------
    (anode_ids, cathode_ids)
        Lists of electrode IDs with significant signals.
    """
    an_ids = list(result.an_collector_id)
    ca_ids = list(result.ca_collector_id)

    # Compute peak absolute signal for each electrode
    an_peaks = [np.max(np.abs(wf)) if len(wf) > 0 else 0.0
                for wf in result.an_vs_time]
    ca_peaks = [np.max(np.abs(wf)) if len(wf) > 0 else 0.0
                for wf in result.ca_vs_time]

    global_max = max(
        max(an_peaks, default=0.0),
        max(ca_peaks, default=0.0),
    )

    if global_max == 0.0:
        return an_ids, ca_ids

    cutoff = threshold * global_max

    selected_anodes = [an_ids[i] for i, p in enumerate(an_peaks) if p >= cutoff]
    selected_cathodes = [ca_ids[i] for i, p in enumerate(ca_peaks) if p >= cutoff]

    return selected_anodes, selected_cathodes

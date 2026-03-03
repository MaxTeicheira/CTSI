"""CTSI Toolkit — unified interface for CZT detector simulation."""

from .config import DetectorSpec, CtsiConfig, load_detector_spec, load_ctsi_config
from .parser import SimulationResult, parse_interactive_output

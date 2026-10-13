"""Molecular dynamics: constant-temperature and heating-ramp protocols.

Stage contract
--------------
Inputs:  minimised holo structure
Outputs: trajectories, rmsf.tsv, timeseries.tsv, qc.json
Config:  dynamics.*
"""

from __future__ import annotations

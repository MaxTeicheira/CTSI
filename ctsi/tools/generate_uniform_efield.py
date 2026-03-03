#!/usr/bin/env python3
"""
Generate a simple uniform electric field for testing CTSI.
E = -600V / 5mm = -120 V/mm = -120,000 V/m
"""

import numpy as np

# Grid parameters (matching ctsi.cfg)
X_MIN, X_MAX = -20.0, 20.0  # mm
Z_MIN, Z_MAX = 0.0, 5.0     # mm
GRID_SPACING = 0.01         # mm

NX = int((X_MAX - X_MIN) / GRID_SPACING) + 1
NZ = int((Z_MAX - Z_MIN) / GRID_SPACING) + 1

x = np.linspace(X_MIN, X_MAX, NX)
z = np.linspace(Z_MIN, Z_MAX, NZ)

# Uniform field: Ez = -120,000 V/m, Ex = 0
# This corresponds to -600V across 5mm
Ez_uniform = -120000.0  # V/m
Ex_uniform = 0.0

print(f"Grid: {NX} x {NZ} = {NX * NZ:,} points")
print(f"Ez = {Ez_uniform} V/m = {Ez_uniform/100} V/cm")

output_file = "./efield_uniform_test.txt"
print(f"Writing to {output_file}...")

with open(output_file, 'w') as f:
    f.write(f"{NX} {NZ} 1 0\n")
    for i in range(NX):
        for j in range(NZ):
            f.write(f"{x[i]:.6g} {z[j]:.6g} {Ex_uniform:.6g} {Ez_uniform:.6g}\n")

print("Done!")

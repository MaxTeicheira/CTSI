#!/usr/bin/env python3
"""Visualize the v3 e-field solution."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load the e-field file
print("Loading e-field data...")
data = []
with open("efield_L5mm_W40mm_600V.txt", 'r') as f:
    header = f.readline()
    nx, nz = int(header.split()[0]), int(header.split()[1])
    for line in f:
        parts = line.split()
        if len(parts) >= 4:
            data.append([float(p) for p in parts[:4]])

data = np.array(data)
x = data[:, 0].reshape(nx, nz)[:, 0]
z = data[:, 1].reshape(nx, nz)[0, :]
Ex = data[:, 2].reshape(nx, nz)
Ez = data[:, 3].reshape(nx, nz)

# Compute potential by integration (Ez = -dV/dz)
dz = z[1] - z[0]
V = np.zeros_like(Ez)
V[:, 0] = -600.0  # Cathode voltage
for j in range(1, len(z)):
    V[:, j] = V[:, j-1] - Ez[:, j-1] * dz / 1000  # Convert V/m to V/mm

print(f"Grid: {nx} x {nz}")
print(f"x range: {x.min():.1f} to {x.max():.1f} mm")
print(f"z range: {z.min():.1f} to {z.max():.1f} mm")

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Potential colormap (zoomed to center)
ax1 = axes[0, 0]
x_mask = (x >= -5) & (x <= 5)
x_zoom = x[x_mask]
V_zoom = V[x_mask, :]

im1 = ax1.pcolormesh(x_zoom, z, V_zoom.T, shading='auto', cmap='RdBu_r')
ax1.set_xlabel('x (mm)')
ax1.set_ylabel('z (mm)')
ax1.set_title('Potential V (zoomed: -5 to 5 mm)')
plt.colorbar(im1, ax=ax1, label='V (volts)')

# Mark electrode positions
anode_centers = np.arange(-19, 20, 1.0)
for ac in anode_centers:
    if -5 <= ac <= 5:
        ax1.axvline(ac, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
ax1.axhline(5.0, color='g', linewidth=2, label='Anode plane')
ax1.axhline(0.0, color='r', linewidth=2, label='Cathode plane')

# 2. Ez field (zoomed)
ax2 = axes[0, 1]
Ez_zoom = Ez[x_mask, :] / 100  # Convert to V/cm
im2 = ax2.pcolormesh(x_zoom, z, Ez_zoom.T, shading='auto', cmap='viridis',
                      vmin=-2000, vmax=0)
ax2.set_xlabel('x (mm)')
ax2.set_ylabel('z (mm)')
ax2.set_title('Ez field (zoomed: -5 to 5 mm)')
plt.colorbar(im2, ax=ax2, label='Ez (V/cm)')

# 3. Potential at anode plane (z=L)
ax3 = axes[1, 0]
ax3.plot(x, V[:, -1], 'b-', linewidth=1)
ax3.set_xlabel('x (mm)')
ax3.set_ylabel('V (volts)')
ax3.set_title('Potential at Anode Plane (z = 5mm)')
ax3.axhline(0, color='k', linestyle='--', alpha=0.5)
ax3.set_xlim(-10, 10)

# Mark anode positions
for ac in anode_centers:
    if -10 <= ac <= 10:
        ax3.axvline(ac, color='g', alpha=0.3, linewidth=0.5)

ax3.text(0.02, 0.98, f'On anode (x=2mm): V = {V[np.argmin(np.abs(x-2)), -1]:.1f}V\n'
                      f'In gap (x=0.5mm): V = {V[np.argmin(np.abs(x-0.5)), -1]:.1f}V',
         transform=ax3.transAxes, verticalalignment='top', fontsize=10,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 4. Ez profiles at different x positions
ax4 = axes[1, 1]
positions = [0.0, 0.5, 1.0, 2.0]
colors = ['r', 'orange', 'g', 'b']
for pos, color in zip(positions, colors):
    idx = np.argmin(np.abs(x - pos))
    ax4.plot(z, Ez[idx, :]/100, color=color, linewidth=1.5,
             label=f'x={pos}mm')

ax4.set_xlabel('z (mm)')
ax4.set_ylabel('Ez (V/cm)')
ax4.set_title('Ez vs z at different x positions')
ax4.axhline(-1200, color='k', linestyle='--', alpha=0.5, label='Uniform: -1200 V/cm')
ax4.legend(loc='lower left')
ax4.set_ylim(-3000, 0)

plt.tight_layout()
plt.savefig('efield_v3_analysis.png', dpi=150, bbox_inches='tight')
print("\nSaved: efield_v3_analysis.png")

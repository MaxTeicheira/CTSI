#!/usr/bin/env python3
"""
Generate 3D weighting potential files for CTSI simulation.

Weighting potential φ_w satisfies Laplace equation with:
- φ_w = 1 on the electrode of interest
- φ_w = 0 on all other electrodes and boundaries

For cross-strip detector:
- Anode weighting potential: solve with single anode pixel = 1, all else = 0
- Cathode weighting potential: solve with single cathode strip = 1, all else = 0
"""

import numpy as np
from scipy.ndimage import laplace
import time

# Detector geometry (matching detectorSpec.txt)
L = 0.5          # Thickness in cm
W = 4.0          # Width in cm (40mm)

NUM_ANODES = 39
ANODE_PITCH = 0.1    # cm (1mm)
ANODE_WIDTH = 0.01   # cm (100um)

NUM_CATHODES = 8
CATHODE_PITCH = 0.5  # cm (5mm)
CATHODE_WIDTH = 0.49 # cm (4900um)

def solve_laplace_3d(phi, fixed_mask, tol=1e-5, max_iter=10000):
    """
    Solve Laplace equation using Jacobi iteration.

    Args:
        phi: Initial potential array (3D)
        fixed_mask: Boolean array, True where potential is fixed
        tol: Convergence tolerance
        max_iter: Maximum iterations

    Returns:
        Converged potential array
    """
    nx, ny, nz = phi.shape

    for iteration in range(max_iter):
        phi_old = phi.copy()

        # Jacobi update for interior points
        phi[1:-1, 1:-1, 1:-1] = (
            phi_old[2:, 1:-1, 1:-1] + phi_old[:-2, 1:-1, 1:-1] +
            phi_old[1:-1, 2:, 1:-1] + phi_old[1:-1, :-2, 1:-1] +
            phi_old[1:-1, 1:-1, 2:] + phi_old[1:-1, 1:-1, :-2]
        ) / 6.0

        # Neumann BC on side boundaries (zero normal derivative)
        # Apply BEFORE restoring fixed values so fixed BCs take precedence
        phi[0, :, :] = phi[1, :, :]
        phi[-1, :, :] = phi[-2, :, :]
        phi[:, 0, :] = phi[:, 1, :]
        phi[:, -1, :] = phi[:, -2, :]

        # Restore fixed values (AFTER Neumann so Dirichlet takes priority)
        phi[fixed_mask] = phi_old[fixed_mask]

        # Check convergence
        max_change = np.max(np.abs(phi - phi_old))

        if iteration % 500 == 0:
            print(f"  Iteration {iteration}: max_change = {max_change:.2e}")

        if max_change < tol:
            print(f"  Converged after {iteration} iterations")
            return phi

    print(f"  Warning: Did not converge after {max_iter} iterations")
    return phi


def generate_anode_weighting_potential(output_file, grid_spacing=0.02):
    """
    Generate weighting potential for central anode pixel.

    The anode pixel is at x=0, running along y direction.
    φ_w = 1 on this pixel at z=L, φ_w = 0 elsewhere.

    IMPORTANT: CTSI uses abs(x), abs(y) for lookup, so we only need
    to generate the positive quadrant (x >= 0, y >= 0).
    """
    print("\n" + "="*60)
    print("Generating ANODE weighting potential")
    print("="*60)

    # Grid setup - only positive quadrant (CTSI uses symmetry)
    # Need to cover at least half anode pitch (0.5mm) with margin
    # Also need enough y range to avoid edge effects
    x_max = 0.1   # cm (1mm from center, covers half anode pitch)
    y_max = 0.5   # cm (5mm along anode, to handle any y position)

    nx = int(x_max / grid_spacing) + 1
    ny = int(y_max / grid_spacing) + 1
    nz = int(L / grid_spacing) + 1

    x = np.linspace(0, x_max, nx)
    y = np.linspace(0, y_max, ny)
    z = np.linspace(0, L, nz)

    print(f"Grid: {nx} x {ny} x {nz} = {nx*ny*nz:,} points")
    print(f"x range: {x[0]:.3f} to {x[-1]:.3f} cm (positive only, CTSI uses symmetry)")
    print(f"y range: {y[0]:.3f} to {y[-1]:.3f} cm (positive only, CTSI uses symmetry)")
    print(f"z range: {z[0]:.3f} to {z[-1]:.3f} cm")

    # Initialize potential
    phi = np.zeros((nx, ny, nz))
    fixed_mask = np.zeros((nx, ny, nz), dtype=bool)

    # NOTE: CTSI uses z_lookup = L - zPos for anode. But empirically, the OLD orientation
    # (φ_w=1 at z=L in file) produces signals while the "correct" orientation doesn't.
    # This might be due to how the signal formula works with qMobile×φ_w sums.

    # Boundary conditions at z=L (anode plane)
    # Central anode pixel: x <= ANODE_WIDTH/2, all y (since x >= 0)
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            if xi <= ANODE_WIDTH / 2:
                # On the central anode pixel
                phi[i, j, -1] = 1.0
                fixed_mask[i, j, -1] = True
            else:
                # On gaps between anodes (at anode plane)
                # For weighting potential, gaps are at 0V
                phi[i, j, -1] = 0.0
                fixed_mask[i, j, -1] = True

    # Boundary conditions at z=0 (cathode plane)
    # All cathodes at φ_w = 0
    phi[:, :, 0] = 0.0
    fixed_mask[:, :, 0] = True

    # Neumann BC at x=0 plane (symmetry: dφ/dx = 0)
    # This is automatically handled by the solver's side boundary treatment

    print("\nSolving Laplace equation...")
    t0 = time.time()
    phi = solve_laplace_3d(phi, fixed_mask, tol=1e-5, max_iter=20000)
    print(f"Solve time: {time.time() - t0:.1f} s")

    # Write output file
    print(f"\nWriting to {output_file}...")
    with open(output_file, 'w') as f:
        # Header: dimensions
        f.write(f"{nx} {ny} {nz}\n")

        # X coordinates
        for xi in x:
            f.write(f"{xi:.6g} ")
        f.write("\n")

        # Y coordinates
        for yi in y:
            f.write(f"{yi:.6g} ")
        f.write("\n")

        # Z coordinates
        for zi in z:
            f.write(f"{zi:.6g} ")
        f.write("\n")

        # Potential values (x-major, y, z innermost)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    f.write(f"{phi[i,j,k]:.6g}\n")

    print(f"Done! φ_w range: {phi.min():.4f} to {phi.max():.4f}")
    return phi, x, y, z


def generate_cathode_weighting_potential(output_file, grid_spacing=0.02):
    """
    Generate weighting potential for central cathode strip.

    The cathode strip runs along detector x direction, with width in y direction.
    φ_w = 1 on this strip at z=0, φ_w = 0 elsewhere.

    IMPORTANT: CTSI coordinate mapping for cathode weighting potential:
    - file x-axis = detector y direction (perpendicular to strip, across width)
    - file y-axis = detector x direction (along strip length)
    - file z-axis = detector z direction (depth)

    CTSI uses abs(x), abs(y) for lookup, so we only need positive quadrant.
    """
    print("\n" + "="*60)
    print("Generating CATHODE weighting potential")
    print("="*60)

    # Grid setup - only positive quadrant (CTSI uses symmetry)
    # File x = detector y (across cathode width)
    # File y = detector x (along cathode length)
    # Cathode pitch = 5mm = 0.5cm, so need to cover at least half (0.25cm)
    x_max = 0.3   # cm (3mm, covers half cathode pitch with margin)
    y_max = 0.1   # cm (1mm, cathode runs along x so less variation needed)

    nx = int(x_max / grid_spacing) + 1
    ny = int(y_max / grid_spacing) + 1
    nz = int(L / grid_spacing) + 1

    x = np.linspace(0, x_max, nx)
    y = np.linspace(0, y_max, ny)
    z = np.linspace(0, L, nz)

    print(f"Grid: {nx} x {ny} x {nz} = {nx*ny*nz:,} points")
    print(f"x range: {x[0]:.3f} to {x[-1]:.3f} cm (= detector y, across strip width)")
    print(f"y range: {y[0]:.3f} to {y[-1]:.3f} cm (= detector x, along strip)")
    print(f"z range: {z[0]:.3f} to {z[-1]:.3f} cm")

    # Initialize potential
    phi = np.zeros((nx, ny, nz))
    fixed_mask = np.zeros((nx, ny, nz), dtype=bool)

    # Boundary conditions at z=0 (cathode plane)
    # Central cathode strip: file_x <= CATHODE_WIDTH/2 (since file_x >= 0)
    # This corresponds to |detector_y| <= CATHODE_WIDTH/2
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            if xi <= CATHODE_WIDTH / 2:
                # On the central cathode strip
                phi[i, j, 0] = 1.0
                fixed_mask[i, j, 0] = True
            else:
                # On gaps between cathodes
                phi[i, j, 0] = 0.0
                fixed_mask[i, j, 0] = True

    # Boundary conditions at z=L (anode plane)
    # All anodes at φ_w = 0
    phi[:, :, -1] = 0.0
    fixed_mask[:, :, -1] = True

    print("\nSolving Laplace equation...")
    t0 = time.time()
    phi = solve_laplace_3d(phi, fixed_mask, tol=1e-5, max_iter=20000)
    print(f"Solve time: {time.time() - t0:.1f} s")

    # Write output file
    print(f"\nWriting to {output_file}...")
    with open(output_file, 'w') as f:
        # Header: dimensions
        f.write(f"{nx} {ny} {nz}\n")

        # X coordinates
        for xi in x:
            f.write(f"{xi:.6g} ")
        f.write("\n")

        # Y coordinates
        for yi in y:
            f.write(f"{yi:.6g} ")
        f.write("\n")

        # Z coordinates
        for zi in z:
            f.write(f"{zi:.6g} ")
        f.write("\n")

        # Potential values
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    f.write(f"{phi[i,j,k]:.6g}\n")

    print(f"Done! φ_w range: {phi.min():.4f} to {phi.max():.4f}")
    return phi, x, y, z


if __name__ == "__main__":
    print("="*60)
    print("WEIGHTING POTENTIAL GENERATOR")
    print("="*60)
    print(f"\nDetector: {W*10:.0f}mm x {L*10:.0f}mm")
    print(f"Anodes: {NUM_ANODES} pixels, {ANODE_WIDTH*10000:.0f}µm wide, {ANODE_PITCH*10:.0f}mm pitch")
    print(f"Cathodes: {NUM_CATHODES} strips, {CATHODE_WIDTH*10000:.0f}µm wide, {CATHODE_PITCH*10:.0f}mm pitch")

    # Generate anode weighting potential
    # Use finer grid (0.002 cm = 20 µm) to resolve 100 µm anode width
    generate_anode_weighting_potential("phi_anode.txt", grid_spacing=0.002)

    # Generate cathode weighting potential
    # Can use coarser grid since cathode is wide (4900 µm)
    generate_cathode_weighting_potential("phi_cathode.txt", grid_spacing=0.01)

    print("\n" + "="*60)
    print("Complete! Generated files:")
    print("  - phi_anode.txt")
    print("  - phi_cathode.txt")
    print("="*60)

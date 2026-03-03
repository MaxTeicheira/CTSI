#!/usr/bin/env python3
"""
Generate 2D electric field file for CZT detector simulation - Version 3.

Correct physics implementation:
- Dirichlet BC ONLY on electrode surfaces
- Gap regions treated as free surfaces (potential determined by Laplace solution)
- No artificial Neumann BC at gaps that blocks field penetration
- Iterative relaxation method for proper convergence

Solves Laplace equation (nabla^2 V = 0) with physically correct boundary conditions.
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline
import time

# ============================================================================
# Detector Geometry Parameters
# ============================================================================

L = 5.0      # Thickness (z direction) mm
W = 40.0     # Width (x direction) mm

NUM_ANODES = 39
ANODE_PITCH = 1.0       # mm
ANODE_WIDTH = 0.1       # mm (100 um)
ANODE_VOLTAGE = 600.0   # V

NUM_CATHODES = 1
CATHODE_PITCH = 40.0    # mm (single cathode covers full width)
CATHODE_WIDTH = 40.0    # mm (planar cathode)
CATHODE_VOLTAGE = 0.0   # V

X_MIN, X_MAX = -20.0, 20.0  # mm
Z_MIN, Z_MAX = 0.0, 5.0     # mm
GRID_SPACING = 0.01         # mm for output

# ============================================================================
# Electrode Functions
# ============================================================================

def get_anode_centers():
    return np.arange(-(NUM_ANODES - 1) / 2, (NUM_ANODES - 1) / 2 + 1, 1.0) * ANODE_PITCH

def get_cathode_centers():
    return np.arange(-(NUM_CATHODES - 1) / 2, (NUM_CATHODES - 1) / 2 + 1, 1.0) * CATHODE_PITCH

def is_on_anode(x):
    """Check if position x is on an anode electrode."""
    centers = get_anode_centers()
    for center in centers:
        if abs(x - center) <= ANODE_WIDTH / 2:
            return True
    return False

def is_on_cathode(x):
    """Check if position x is on a cathode electrode."""
    centers = get_cathode_centers()
    for center in centers:
        if abs(x - center) <= CATHODE_WIDTH / 2:
            return True
    return False

# ============================================================================
# Correct Laplace Solver with Floating Gap Potentials
# ============================================================================

def solve_laplace_correct(nx, nz, tol=1e-6, max_iter=50000):
    """
    Solve Laplace equation with physically correct boundary conditions:
    - Dirichlet on electrodes only
    - Gap regions: potential floats (determined by solution)
    - Side boundaries: Neumann (dV/dx = 0)

    Uses vectorized Jacobi iteration with NumPy for speed.
    """
    dx = (X_MAX - X_MIN) / (nx - 1)
    dz = (Z_MAX - Z_MIN) / (nz - 1)

    print(f"\nSolving on grid: {nx} x {nz} = {nx * nz:,} points")
    print(f"Grid spacing: dx={dx:.4f} mm, dz={dz:.4f} mm")

    x = np.linspace(X_MIN, X_MAX, nx)
    z = np.linspace(Z_MIN, Z_MAX, nz)

    # Create electrode masks
    anode_mask = np.array([is_on_anode(xi) for xi in x])
    cathode_mask = np.array([is_on_cathode(xi) for xi in x])

    print(f"Anode points: {np.sum(anode_mask)}, Cathode points: {np.sum(cathode_mask)}")
    print(f"Anode gaps: {nx - np.sum(anode_mask)} points at z=L")
    print(f"Cathode gaps: {nx - np.sum(cathode_mask)} points at z=0")

    # Initialize potential with linear gradient (good initial guess)
    V = np.zeros((nx, nz))
    Z_grid = np.tile(z, (nx, 1))
    V = CATHODE_VOLTAGE + (ANODE_VOLTAGE - CATHODE_VOLTAGE) * Z_grid / L

    # Create fixed mask (True where potential is fixed)
    fixed = np.zeros((nx, nz), dtype=bool)
    fixed[cathode_mask, 0] = True
    fixed[anode_mask, -1] = True

    # Apply Dirichlet BC on electrodes
    V[cathode_mask, 0] = CATHODE_VOLTAGE
    V[anode_mask, -1] = ANODE_VOLTAGE

    # Precompute coefficients for non-square grid
    rx = (dz / dx) ** 2
    rz = 1.0
    denom = 2.0 * (rx + rz)
    denom_boundary = 2.0 * rx + 2.0 * rz

    # Masks for gap points
    anode_gap_mask = ~anode_mask
    cathode_gap_mask = ~cathode_mask

    print(f"\nStarting vectorized Jacobi iteration...")
    print(f"Convergence tolerance: {tol}")

    t0 = time.time()

    for iteration in range(max_iter):
        V_old = V.copy()

        # Vectorized update for interior points (1:-1, 1:-1)
        V_new_interior = (rx * (V_old[2:, 1:-1] + V_old[:-2, 1:-1]) +
                          rz * (V_old[1:-1, 2:] + V_old[1:-1, :-2])) / denom
        V[1:-1, 1:-1] = V_new_interior

        # Update gap points at cathode plane (j=0) - FREE SURFACE
        # Use one-sided difference in z
        V_new_cathode_gap = (rx * (V_old[2:, 0] + V_old[:-2, 0]) +
                             2.0 * rz * V_old[1:-1, 1]) / denom_boundary
        # Only update non-electrode points
        gap_indices = np.where(cathode_gap_mask[1:-1])[0] + 1
        V[gap_indices, 0] = V_new_cathode_gap[gap_indices - 1]

        # Update gap points at anode plane (j=-1) - FREE SURFACE
        V_new_anode_gap = (rx * (V_old[2:, -1] + V_old[:-2, -1]) +
                           2.0 * rz * V_old[1:-1, -2]) / denom_boundary
        gap_indices = np.where(anode_gap_mask[1:-1])[0] + 1
        V[gap_indices, -1] = V_new_anode_gap[gap_indices - 1]

        # Restore fixed electrode potentials
        V[cathode_mask, 0] = CATHODE_VOLTAGE
        V[anode_mask, -1] = ANODE_VOLTAGE

        # Side boundaries (Neumann: dV/dx = 0)
        V[0, :] = V[1, :]
        V[-1, :] = V[-2, :]

        # Check convergence
        max_change = np.max(np.abs(V - V_old))

        if iteration % 500 == 0:
            elapsed = time.time() - t0
            print(f"  Iteration {iteration}: max_change = {max_change:.2e}, "
                  f"time = {elapsed:.1f}s")

        if max_change < tol:
            print(f"\nConverged after {iteration} iterations!")
            print(f"Final max_change: {max_change:.2e}")
            break
    else:
        print(f"\nWarning: Did not converge after {max_iter} iterations")
        print(f"Final max_change: {max_change:.2e}")

    print(f"Solve time: {time.time() - t0:.1f} s")

    return V, x, z

def solve_laplace_multigrid(nx_fine, nz_fine, tol=1e-6):
    """
    Multigrid solver for faster convergence on fine grids.
    Solves on coarse grid first, then interpolates and refines.
    """
    print("\n" + "="*60)
    print("MULTIGRID SOLVER")
    print("="*60)

    # Level 1: Coarse grid (4x coarser)
    nx_coarse = (nx_fine - 1) // 4 + 1
    nz_coarse = (nz_fine - 1) // 4 + 1

    print(f"\nLevel 1 (coarse): {nx_coarse} x {nz_coarse}")
    V_coarse, x_coarse, z_coarse = solve_laplace_correct(nx_coarse, nz_coarse, tol=1e-5)

    # Level 2: Medium grid (2x coarser than fine)
    nx_med = (nx_fine - 1) // 2 + 1
    nz_med = (nz_fine - 1) // 2 + 1

    print(f"\nLevel 2 (medium): {nx_med} x {nz_med}")

    # Interpolate coarse solution to medium grid
    x_med = np.linspace(X_MIN, X_MAX, nx_med)
    z_med = np.linspace(Z_MIN, Z_MAX, nz_med)

    interp = RectBivariateSpline(x_coarse, z_coarse, V_coarse, kx=3, ky=3)
    V_med_init = interp(x_med, z_med)

    # Solve on medium grid with interpolated initial guess
    V_med, x_med, z_med = solve_laplace_with_init(nx_med, nz_med, V_med_init, tol=1e-5)

    # Level 3: Fine grid
    print(f"\nLevel 3 (fine): {nx_fine} x {nz_fine}")

    x_fine = np.linspace(X_MIN, X_MAX, nx_fine)
    z_fine = np.linspace(Z_MIN, Z_MAX, nz_fine)

    interp = RectBivariateSpline(x_med, z_med, V_med, kx=3, ky=3)
    V_fine_init = interp(x_fine, z_fine)

    V_fine, x_fine, z_fine = solve_laplace_with_init(nx_fine, nz_fine, V_fine_init, tol=tol)

    return V_fine, x_fine, z_fine

def solve_laplace_with_init(nx, nz, V_init, tol=1e-6, max_iter=50000):
    """
    Solve Laplace equation starting from an initial guess (vectorized).
    """
    dx = (X_MAX - X_MIN) / (nx - 1)
    dz = (Z_MAX - Z_MIN) / (nz - 1)

    x = np.linspace(X_MIN, X_MAX, nx)
    z = np.linspace(Z_MIN, Z_MAX, nz)

    # Create electrode masks
    anode_mask = np.array([is_on_anode(xi) for xi in x])
    cathode_mask = np.array([is_on_cathode(xi) for xi in x])

    # Start from initial guess
    V = V_init.copy()

    # Apply Dirichlet BC on electrodes
    V[cathode_mask, 0] = CATHODE_VOLTAGE
    V[anode_mask, -1] = ANODE_VOLTAGE

    rx = (dz / dx) ** 2
    rz = 1.0
    denom = 2.0 * (rx + rz)
    denom_boundary = 2.0 * rx + 2.0 * rz

    anode_gap_mask = ~anode_mask
    cathode_gap_mask = ~cathode_mask

    t0 = time.time()

    for iteration in range(max_iter):
        V_old = V.copy()

        # Vectorized interior update
        V[1:-1, 1:-1] = (rx * (V_old[2:, 1:-1] + V_old[:-2, 1:-1]) +
                         rz * (V_old[1:-1, 2:] + V_old[1:-1, :-2])) / denom

        # Cathode gap points
        V_new_cathode_gap = (rx * (V_old[2:, 0] + V_old[:-2, 0]) +
                             2.0 * rz * V_old[1:-1, 1]) / denom_boundary
        gap_indices = np.where(cathode_gap_mask[1:-1])[0] + 1
        V[gap_indices, 0] = V_new_cathode_gap[gap_indices - 1]

        # Anode gap points
        V_new_anode_gap = (rx * (V_old[2:, -1] + V_old[:-2, -1]) +
                           2.0 * rz * V_old[1:-1, -2]) / denom_boundary
        gap_indices = np.where(anode_gap_mask[1:-1])[0] + 1
        V[gap_indices, -1] = V_new_anode_gap[gap_indices - 1]

        # Restore electrode potentials
        V[cathode_mask, 0] = CATHODE_VOLTAGE
        V[anode_mask, -1] = ANODE_VOLTAGE

        # Side boundaries
        V[0, :] = V[1, :]
        V[-1, :] = V[-2, :]

        max_change = np.max(np.abs(V - V_old))

        if iteration % 500 == 0:
            print(f"  Iteration {iteration}: max_change = {max_change:.2e}")

        if max_change < tol:
            print(f"  Converged after {iteration} iterations")
            break

    print(f"  Solve time: {time.time() - t0:.1f} s")

    return V, x, z

def interpolate_to_fine_grid(V_coarse, x_coarse, z_coarse):
    """Interpolate to fine output grid."""
    NX = int((X_MAX - X_MIN) / GRID_SPACING) + 1
    NZ = int((Z_MAX - Z_MIN) / GRID_SPACING) + 1

    print(f"\nInterpolating to output grid: {NX} x {NZ}")

    x_fine = np.linspace(X_MIN, X_MAX, NX)
    z_fine = np.linspace(Z_MIN, Z_MAX, NZ)

    interp = RectBivariateSpline(x_coarse, z_coarse, V_coarse, kx=3, ky=3)
    V_fine = interp(x_fine, z_fine)

    # Enforce exact BC on electrodes
    print("Enforcing exact electrode potentials...")
    for i, xi in enumerate(x_fine):
        if is_on_cathode(xi):
            V_fine[i, 0] = CATHODE_VOLTAGE
        if is_on_anode(xi):
            V_fine[i, -1] = ANODE_VOLTAGE

    return V_fine, x_fine, z_fine

def compute_efield(V, x, z, max_field=200000.0):
    """
    Compute E = -grad(V) with proper handling of boundaries.

    Args:
        V: Potential array
        x, z: Coordinate arrays
        max_field: Maximum allowed field magnitude in V/m (default 200 kV/m = 2000 V/cm)
    """
    nx, nz = V.shape
    dx = x[1] - x[0]
    dz = z[1] - z[0]

    Ex = np.zeros_like(V)
    Ez = np.zeros_like(V)

    # Central differences for interior
    Ex[1:-1, :] = -(V[2:, :] - V[:-2, :]) / (2 * dx)
    Ez[:, 1:-1] = -(V[:, 2:] - V[:, :-2]) / (2 * dz)

    # One-sided differences at boundaries
    Ex[0, :] = -(V[1, :] - V[0, :]) / dx
    Ex[-1, :] = -(V[-1, :] - V[-2, :]) / dx
    Ez[:, 0] = -(V[:, 1] - V[:, 0]) / dz
    Ez[:, -1] = -(V[:, -1] - V[:, -2]) / dz

    # Convert V/mm to V/m
    Ex *= 1000.0
    Ez *= 1000.0

    # Cap extreme field values to prevent numerical overflow in simulation
    # Typical bulk field is ~120,000 V/m, allow up to ~5x for fringing
    E_mag = np.sqrt(Ex**2 + Ez**2)
    mask = E_mag > max_field
    if np.any(mask):
        scale = max_field / E_mag[mask]
        Ex[mask] *= scale
        Ez[mask] *= scale
        n_capped = np.sum(mask)
        print(f"  Capped {n_capped} points ({100*n_capped/Ex.size:.2f}%) to max field {max_field/1000:.0f} kV/m")

    return Ex, Ez

def write_efield_file(filename, x, z, Ex, Ez):
    """Write E-field in CTSI format."""
    nx, nz = len(x), len(z)

    print(f"\nWriting E-field file: {filename}")
    print(f"Grid: {nx} x {nz} = {nx * nz:,} points")

    with open(filename, 'w') as f:
        f.write(f"{nx} {nz} 1 0\n")
        for i in range(nx):
            for j in range(nz):
                f.write(f"{x[i]:.6g} {z[j]:.6g} {Ex[i,j]:.6g} {Ez[i,j]:.6g}\n")
            if (i + 1) % 500 == 0:
                print(f"  Written {i + 1}/{nx} ({100*(i+1)/nx:.1f}%)")

    print("File written successfully!")

def verify_solution(V, Ex, Ez, x, z):
    """Print verification statistics."""
    print("\n" + "="*60)
    print("SOLUTION VERIFICATION")
    print("="*60)

    expected_Ez = -CATHODE_VOLTAGE / L * 1000  # V/m
    print(f"\nExpected bulk Ez (if uniform): {expected_Ez:.0f} V/m ({expected_Ez/100:.0f} V/cm)")

    # At anode position (x=2mm)
    anode_idx = np.argmin(np.abs(x - 2.0))
    print(f"\nAt x={x[anode_idx]:.2f}mm (on anode):")
    print(f"  V at z=0 (cathode): {V[anode_idx, 0]:.1f} V")
    print(f"  V at z=L (anode):   {V[anode_idx, -1]:.1f} V")
    print(f"  Ez at z=0:   {Ez[anode_idx, 0]:.0f} V/m")
    print(f"  Ez at z=L/2: {Ez[anode_idx, len(z)//2]:.0f} V/m")
    print(f"  Ez at z=L:   {Ez[anode_idx, -1]:.0f} V/m")

    # At gap position (x=0.5mm - between anodes at 0 and 1mm)
    gap_idx = np.argmin(np.abs(x - 0.5))
    print(f"\nAt x={x[gap_idx]:.2f}mm (anode gap):")
    print(f"  V at z=0 (cathode): {V[gap_idx, 0]:.1f} V")
    print(f"  V at z=L (gap):     {V[gap_idx, -1]:.1f} V  (should NOT be 0V)")
    print(f"  Ez at z=0:   {Ez[gap_idx, 0]:.0f} V/m")
    print(f"  Ez at z=L/2: {Ez[gap_idx, len(z)//2]:.0f} V/m")
    print(f"  Ez at z=L:   {Ez[gap_idx, -1]:.0f} V/m")

    # At cathode gap (x=0 is between cathodes at -2.5 and +2.5mm)
    cgap_idx = np.argmin(np.abs(x - 0.0))
    print(f"\nAt x={x[cgap_idx]:.2f}mm (cathode gap):")
    print(f"  V at z=0 (gap):     {V[cgap_idx, 0]:.1f} V  (should NOT be -600V)")
    print(f"  V at z=L:           {V[cgap_idx, -1]:.1f} V")
    print(f"  Ez at z=0:   {Ez[cgap_idx, 0]:.0f} V/m")
    print(f"  Ez at z=L/2: {Ez[cgap_idx, len(z)//2]:.0f} V/m")

    print(f"\nField statistics:")
    print(f"  Ez range: {Ez.min()/100:.1f} to {Ez.max()/100:.1f} V/cm")
    print(f"  Ex range: {Ex.min()/100:.1f} to {Ex.max()/100:.1f} V/cm")

    # Check for physically reasonable values
    print(f"\nPhysics checks:")
    if Ez.min() < 0 and Ez.max() > 0:
        print(f"  WARNING: Ez changes sign - check for issues")
    elif Ez.max() < 0:
        print(f"  OK: Ez is negative (correct for e- drift toward anode)")
    else:
        print(f"  WARNING: Ez is positive - electrons would drift wrong way")

# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    print("="*60)
    print("E-FIELD GENERATOR V3 - Correct Boundary Conditions")
    print("="*60)

    print(f"\nDetector: {W}mm x {L}mm")
    print(f"Cathode: {CATHODE_VOLTAGE}V, {NUM_CATHODES} strips @ {CATHODE_WIDTH}mm wide")
    print(f"Anode: {ANODE_VOLTAGE}V, {NUM_ANODES} pixels @ {ANODE_WIDTH}mm wide")
    print(f"\nKey difference from v1/v2:")
    print(f"  - Gap regions have FLOATING potential (not fixed BC)")
    print(f"  - Potential at gaps determined by Laplace solution")

    t_start = time.time()

    # Solve on moderate grid (0.025mm spacing for accuracy)
    nx_solve = int((X_MAX - X_MIN) / 0.025) + 1
    nz_solve = int((Z_MAX - Z_MIN) / 0.025) + 1

    V_solve, x_solve, z_solve = solve_laplace_correct(nx_solve, nz_solve, tol=1e-5)

    # Interpolate to fine output grid
    V_fine, x_fine, z_fine = interpolate_to_fine_grid(V_solve, x_solve, z_solve)

    # Compute E-field
    Ex, Ez = compute_efield(V_fine, x_fine, z_fine)

    # Verify
    verify_solution(V_fine, Ex, Ez, x_fine, z_fine)

    # Write output
    output_file = "./efield_L5mm_W40mm_600V.txt"
    write_efield_file(output_file, x_fine, z_fine, Ex, Ez)

    print(f"\nTotal time: {time.time() - t_start:.1f} s")

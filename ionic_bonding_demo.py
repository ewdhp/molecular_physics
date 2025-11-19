"""Ionic bonding demo (simple model)

This script demonstrates the basic theory of ionic bonding using a
pair potential composed of Coulomb attraction and a short-range
Born-Mayer repulsion term:

    V(r) = k * q1 * q2 / r + B * exp(-r / rho)

- Units: SI (meters, Joules) internally; outputs shown in angstroms and eV.
- Parameters chosen for demonstration only (not fitted to a specific crystal).

It finds the equilibrium separation (minimizes V) and plots V(r).
"""

import numpy as np
import matplotlib
# Use non-interactive backend so this runs headless and saves files
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Physical constants
eps0 = 8.8541878128e-12      # vacuum permittivity (F/m)
ke = 1 / (4 * np.pi * eps0)  # Coulomb constant (N·m^2/C^2)
e = 1.602176634e-19         # elementary charge (C)

# Ion charges (for Na+ / Cl- style pair use +1 and -1)
q1 = +1.0 * e
q2 = -1.0 * e

# Born-Mayer repulsion parameters (illustrative)
# V_rep(r) = B * exp(-r / rho)
B = 1e-14      # J (tunable, chosen so repulsion balances Coulomb around a few Angstroms)
rho = 3e-11   # m  (0.3 Angstrom)

# Potential function (SI units)
def V_pair(r):
    """Pair potential (J) for separation r (m)."""
    coulomb = ke * q1 * q2 / r
    rep = B * np.exp(-r / rho)
    return coulomb + rep

# Force (N): -dV/dr
def force_pair(r):
    # analytic derivative of V_pair: d/dr (ke*q1*q2/r) = -ke*q1*q2/r^2
    # d/dr (B exp(-r/rho)) = -B/rho * exp(-r/rho)
    coulomb_deriv = -ke * q1 * q2 / (r ** 2)
    rep_deriv = -B / rho * np.exp(-r / rho)
    # Force on q1 due to q2 is -dV/dr
    return - (coulomb_deriv + rep_deriv)

# r grid (meters) — avoid r=0
r_min = 0.5e-10  # 0.5 Å
r_max = 6.0e-10  # 6.0 Å
r = np.linspace(r_min, r_max, 2000)

V = V_pair(r)
F = force_pair(r)

# Find equilibrium (min potential)
idx_min = np.argmin(V)
r_eq = r[idx_min]
V_eq = V[idx_min]

# Convert units for reporting
m_to_A = 1e10
J_to_eV = 1 / e

print(f"Equilibrium separation r_eq = {r_eq * m_to_A:.3f} Å")
print(f"Potential at r_eq = {V_eq * J_to_eV:.3f} eV")

# Verify force near equilibrium is ~0
F_eq = F[idx_min]
print(f"Force at r_eq = {F_eq:.3e} N (should be ≈ 0)")

# Plot potential (eV) vs r (Å)
plt.figure(figsize=(6,4))
plt.plot(r * m_to_A, V * J_to_eV, label='V(r) (eV)')
plt.axvline(r_eq * m_to_A, color='gray', linestyle='--', alpha=0.6, label=f'r_eq = {r_eq*m_to_A:.3f} Å')
plt.scatter([r_eq * m_to_A], [V_eq * J_to_eV], color='red')
plt.text(r_eq * m_to_A * 1.05, V_eq * J_to_eV, f"  {V_eq*J_to_eV:.2f} eV", va='bottom')
plt.xlabel('Separation r (Å)')
plt.ylabel('Potential energy V (eV)')
plt.title('Simple Ionic Pair Potential (Coulomb + Born–Mayer)')
plt.grid(True)
plt.legend()

out_dir = os.path.dirname(__file__)
out_path = os.path.join(out_dir, 'ionic_bonding_demo.png')
plt.tight_layout()
plt.savefig(out_path, dpi=150)
print(f"Saved potential plot to: {out_path}")

# Save a small data file with r, V, F for further inspection
data_path = os.path.join(out_dir, 'ionic_bonding_demo_data.csv')
np.savetxt(data_path, np.vstack([r * m_to_A, V * J_to_eV, F]).T,
           header='r_A,V_eV,Force_N', delimiter=',')
print(f"Saved data to: {data_path}")

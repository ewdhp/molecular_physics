"""
Quantum Numbers and Electron States Demonstration
==================================================
Demonstrates how quantum numbers (n, l, m_l, m_s) define electron states
in atoms, including orbital shapes, energies, and electron configurations.

Topics covered:
1. Quantum number relationships and constraints
   Introduction: Quantum numbers are discrete values that describe the quantum
   state of an electron in an atom. Unlike classical mechanics where particles
   can have any energy or position, quantum mechanics restricts electrons to
   specific quantum states characterized by four quantum numbers.
   
   Theory: The four quantum numbers arise from solving the Schrödinger equation
   for the hydrogen atom in spherical coordinates:
   - Principal quantum number (n): Determines the energy level and average
     distance from nucleus. From the radial equation: E_n = -13.6 eV/n²
   - Azimuthal/Angular momentum quantum number (l): Determines orbital shape
     and angular momentum L = √(l(l+1))ℏ. Arises from angular part of wave equation.
   - Magnetic quantum number (m_l): Determines orbital orientation in space,
     specifically the z-component of angular momentum: L_z = m_l·ℏ
   - Spin quantum number (m_s): Intrinsic angular momentum of electron (±½),
     discovered through Stern-Gerlach experiment and required by quantum mechanics.

2. Electron configurations for elements
   Introduction: Electron configuration describes how electrons are distributed
   among atomic orbitals in an atom's ground state. This distribution follows
   specific rules that minimize the total energy of the atom.
   
   Theory: The Aufbau principle (German: "building up") states electrons fill
   orbitals in order of increasing energy: 1s → 2s → 2p → 3s → 3p → 4s → 3d...
   This order follows the (n+l) rule: orbitals with lower (n+l) fill first;
   if (n+l) is equal, lower n fills first. This arises from electron-electron
   repulsion and nuclear shielding effects in multi-electron atoms. The Madelung
   rule explains why 4s fills before 3d despite being in a higher shell.

3. Orbital shapes and visualizations
   Introduction: Atomic orbitals are three-dimensional regions where electrons
   are most likely to be found. Each orbital type (s, p, d, f) has a
   characteristic shape determined by the angular part of the wavefunction.
   
   Theory: The complete wavefunction Ψ(r,θ,φ) = R_nl(r)·Y_lm(θ,φ) has two parts:
   - Radial function R_nl(r): Describes probability vs. distance from nucleus,
     containing (n-l-1) radial nodes where probability = 0
   - Angular function Y_lm(θ,φ): Spherical harmonics determining orbital shape,
     containing l angular nodes (nodal planes or cones)
   Shape types: s (spherical), p (dumbbell), d (cloverleaf), f (complex multilobes)

4. Energy level diagrams
   Introduction: Energy levels represent the quantized energies an electron can
   have in an atom. For hydrogen, energy depends only on n; for multi-electron
   atoms, both n and l affect energy due to electron-electron interactions.
   
   Theory: For hydrogen-like atoms (one electron), Bohr model and Schrödinger
   equation both give E_n = -Z²·13.6 eV/n², where Z is nuclear charge. The
   negative energy indicates bound states. As n→∞, E→0 (ionization threshold).
   In multi-electron atoms, orbital energies follow the order: ns < np < nd < nf
   for the same n, due to penetration effect (s electrons penetrate closer to
   nucleus, experiencing less shielding) and electron-electron repulsion.

5. Pauli exclusion principle
   Introduction: The Pauli Exclusion Principle is a fundamental law of quantum
   mechanics stating that no two electrons in an atom can have identical quantum
   numbers. This explains the structure of the periodic table and chemical bonding.
   
   Theory: Formulated by Wolfgang Pauli (1925), this principle arises from the
   antisymmetric nature of fermionic wavefunctions. For fermions (half-integer
   spin particles like electrons): Ψ(r₁,r₂) = -Ψ(r₂,r₁). If two electrons had
   identical quantum numbers, they would occupy the same state, making the
   wavefunction Ψ = -Ψ, thus Ψ = 0 (impossible). Consequence: each orbital
   (defined by n,l,m_l) holds maximum 2 electrons with opposite spins (m_s=±½).

6. Hund's rule demonstration
   Introduction: Hund's rules predict ground state electron configurations for
   atoms, particularly how electrons fill degenerate orbitals (same energy).
   Critical for understanding magnetic properties and chemical reactivity.
   
   Theory: Friedrich Hund's rules (1925) for ground state configurations:
   Rule 1 (Maximum Multiplicity): Electrons occupy degenerate orbitals singly
   with parallel spins before pairing. Reason: parallel spins stay farther apart
   spatially (exchange interaction), reducing electron-electron repulsion energy.
   Rule 2: For given multiplicity, maximum total orbital angular momentum L.
   Rule 3: If subshell ≤ half-filled, minimum J; if > half-filled, maximum J.
   Example: Carbon (2p²) has ↑ ↑ _ configuration (parallel) rather than ↑↓ _ _
   (antiparallel), giving triplet ground state ³P.

Author: ewdhp
Date: November 7, 2025
"""

import numpy as np
import matplotlib
# Try to use an interactive backend
try:
    matplotlib.use('TkAgg')
except:
    try:
        matplotlib.use('Qt5Agg')
    except:
        pass  # Use whatever is available
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm, genlaguerre, factorial
from itertools import product

# Physical constants
a0 = 5.29177e-11  # Bohr radius (m)
hbar = 1.054571817e-34  # Reduced Planck constant (J·s)
me = 9.10938356e-31  # Electron mass (kg)
e = 1.602176634e-19  # Elementary charge (C)
epsilon0 = 8.854187817e-12  # Vacuum permittivity (F/m)


def validate_quantum_numbers(n, l, ml, ms):
    """
    Validate quantum numbers according to quantum mechanical rules
    
    Parameters:
    -----------
    n : int
        Principal quantum number (n ≥ 1)
    l : int
        Azimuthal quantum number (0 ≤ l < n)
    ml : int
        Magnetic quantum number (-l ≤ ml ≤ l)
    ms : float
        Spin quantum number (±1/2)
    
    Returns:
    --------
    bool : True if valid, False otherwise
    """
    if n < 1:
        return False
    if l < 0 or l >= n:
        return False
    if abs(ml) > l:
        return False
    if abs(ms) != 0.5:
        return False
    return True


def orbital_name(n, l):
    """
    Convert quantum numbers to orbital name (e.g., 1s, 2p, 3d)
    
    Parameters:
    -----------
    n : int
        Principal quantum number
    l : int
        Azimuthal quantum number
    
    Returns:
    --------
    str : Orbital name
    """
    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
    return f"{n}{l_names.get(l, '?')}"


def max_electrons_in_orbital(l):
    """
    Maximum number of electrons in an orbital with given l
    
    Parameters:
    -----------
    l : int
        Azimuthal quantum number
    
    Returns:
    --------
    int : Maximum number of electrons (2(2l+1))
    """
    return 2 * (2 * l + 1)


def energy_level_hydrogen(n):
    """
    Energy of electron in hydrogen atom (Bohr model)
    
    Parameters:
    -----------
    n : int
        Principal quantum number
    
    Returns:
    --------
    float : Energy in eV
    """
    # E_n = -13.6 eV / n^2
    return -13.6 / (n**2)


def radial_wavefunction(r, n, l, Z=1):
    """
    Radial part of hydrogen-like wavefunction
    
    Parameters:
    -----------
    r : array
        Radial distance in units of a0 (Bohr radius)
    n : int
        Principal quantum number
    l : int
        Azimuthal quantum number
    Z : int
        Nuclear charge (1 for hydrogen)
    
    Returns:
    --------
    array : Radial wavefunction R_nl(r)
    """
    rho = 2 * Z * r / n
    
    # Normalization constant
    norm = np.sqrt((2*Z/n)**3 * factorial(n-l-1) / (2*n*factorial(n+l)))
    
    # Laguerre polynomial
    laguerre = genlaguerre(n-l-1, 2*l+1)(rho)
    
    # Radial function
    R = norm * np.exp(-rho/2) * rho**l * laguerre
    
    return R


def generate_electron_configuration(Z):
    """
    Generate electron configuration for element with atomic number Z
    
    Parameters:
    -----------
    Z : int
        Atomic number (number of electrons)
    
    Returns:
    --------
    dict : Configuration with orbital occupation
    """
    # Aufbau order (n+l rule)
    orbitals = []
    for n in range(1, 10):
        for l in range(n):
            orbitals.append((n, l))
    
    # Sort by n+l, then by n
    orbitals.sort(key=lambda x: (x[0] + x[1], x[0]))
    
    config = {}
    electrons_left = Z
    
    for n, l in orbitals:
        if electrons_left <= 0:
            break
        
        max_e = max_electrons_in_orbital(l)
        electrons_in_orbital = min(electrons_left, max_e)
        
        orbital_label = orbital_name(n, l)
        config[orbital_label] = electrons_in_orbital
        electrons_left -= electrons_in_orbital
    
    return config


def configuration_string(config):
    """
    Convert configuration dict to standard notation
    
    Parameters:
    -----------
    config : dict
        Electron configuration
    
    Returns:
    --------
    str : Configuration string (e.g., "1s² 2s² 2p⁶")
    """
    superscripts = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    parts = [f"{orb}{count}".translate(superscripts) for orb, count in config.items()]
    return " ".join(parts)


def get_element_name(Z):
    """Get element name from atomic number"""
    elements = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
        9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
        16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 25: "Mn", 26: "Fe"
    }
    return elements.get(Z, f"Z={Z}")


# ============================================================================
# CALCULATIONS
# ============================================================================

print("=" * 80)
print("QUANTUM NUMBERS AND ELECTRON STATES")
print("=" * 80)
print()

print("1. QUANTUM NUMBERS - DEFINITIONS AND CONSTRAINTS")
print("-" * 80)
print()
print("┌─────────┬────────────────────────┬─────────────────┬──────────────────────────┐")
print("│ Symbol  │ Name                   │ Allowed Values  │ Physical Meaning         │")
print("├─────────┼────────────────────────┼─────────────────┼──────────────────────────┤")
print("│ n       │ Principal              │ 1, 2, 3, ...    │ Energy level, size       │")
print("│ l       │ Azimuthal/Angular      │ 0 to (n-1)      │ Orbital shape            │")
print("│ m_l     │ Magnetic               │ -l to +l        │ Orbital orientation      │")
print("│ m_s     │ Spin                   │ +½ or -½        │ Electron spin direction  │")
print("└─────────┴────────────────────────┴─────────────────┴──────────────────────────┘")
print()

print("Orbital designations: l = 0(s), 1(p), 2(d), 3(f), 4(g), 5(h), ...")
print()

# Example quantum number sets
print("2. EXAMPLE QUANTUM NUMBER SETS")
print("-" * 80)
print()

examples = [
    (1, 0, 0, 0.5, "1s orbital, spin up"),
    (2, 0, 0, -0.5, "2s orbital, spin down"),
    (2, 1, -1, 0.5, "2p orbital (m_l=-1), spin up"),
    (2, 1, 0, 0.5, "2p orbital (m_l=0), spin up"),
    (2, 1, 1, 0.5, "2p orbital (m_l=+1), spin up"),
    (3, 2, -2, 0.5, "3d orbital (m_l=-2), spin up"),
]

print(f"{'n':<4} {'l':<4} {'m_l':<6} {'m_s':<6} {'Valid?':<8} {'Description':<40}")
print("-" * 80)
for n, l, ml, ms, description in examples:
    valid = validate_quantum_numbers(n, l, ml, ms)
    print(f"{n:<4} {l:<4} {ml:<6} {ms:<6.1f} {'✓' if valid else '✗':<8} {description:<40}")

print()

# Maximum electrons per orbital type
print("3. MAXIMUM ELECTRONS PER ORBITAL TYPE")
print("-" * 80)
print()
print(f"{'Orbital':<10} {'l value':<10} {'# of orbitals':<15} {'Max electrons':<15}")
print("-" * 80)
for l in range(5):
    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g'}
    num_orbitals = 2*l + 1
    max_e = max_electrons_in_orbital(l)
    print(f"{l_names[l]:<10} {l:<10} {num_orbitals:<15} {max_e:<15}")

print()

# Electron configurations
print("4. ELECTRON CONFIGURATIONS (AUFBAU PRINCIPLE)")
print("-" * 80)
print()

test_elements = [1, 2, 3, 6, 7, 8, 10, 11, 18, 20, 26]

print(f"{'Z':<4} {'Element':<8} {'Configuration':<50}")
print("-" * 80)
for Z in test_elements:
    config = generate_electron_configuration(Z)
    config_str = configuration_string(config)
    element = get_element_name(Z)
    print(f"{Z:<4} {element:<8} {config_str:<50}")

print()

# Detailed analysis for Carbon (Z=6)
print("5. DETAILED ANALYSIS: CARBON (Z = 6)")
print("-" * 80)
print()
print("Carbon has 6 electrons distributed as: 1s² 2s² 2p²")
print()
print("Individual electron quantum numbers:")
print(f"{'Electron':<10} {'n':<4} {'l':<4} {'m_l':<6} {'m_s':<6} {'Orbital':<10} {'Spin':<6}")
print("-" * 80)

carbon_electrons = [
    (1, 1, 0, 0, 0.5, "1s", "↑"),
    (2, 1, 0, 0, -0.5, "1s", "↓"),
    (3, 2, 0, 0, 0.5, "2s", "↑"),
    (4, 2, 0, 0, -0.5, "2s", "↓"),
    (5, 2, 1, -1, 0.5, "2p", "↑"),  # Hund's rule: maximize unpaired spins
    (6, 2, 1, 0, 0.5, "2p", "↑"),   # Second 2p electron in different orbital
]

for e_num, n, l, ml, ms, orb, spin in carbon_electrons:
    print(f"{e_num:<10} {n:<4} {l:<4} {ml:<6} {ms:<6.1f} {orb:<10} {spin:<6}")

print()
print("Note: Electrons 5 and 6 follow Hund's rule - they occupy different")
print("      2p orbitals with parallel spins to minimize electron-electron repulsion.")
print()

# Energy levels for hydrogen
print("6. ENERGY LEVELS (HYDROGEN ATOM)")
print("-" * 80)
print()
print(f"{'n':<4} {'Energy (eV)':<15} {'Orbitals':<30}")
print("-" * 80)
for n in range(1, 6):
    energy = energy_level_hydrogen(n)
    orbitals = [orbital_name(n, l) for l in range(n)]
    print(f"{n:<4} {energy:<15.4f} {', '.join(orbitals):<30}")

print()
print("=" * 80)


# ============================================================================
# VISUALIZATIONS
# ============================================================================

print("\nGenerating visualizations...")
print("Close the plot window to continue.\n")

fig = plt.figure(figsize=(18, 12))
gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

# Color scheme
colors_orbital = {'s': '#FF6B6B', 'p': '#4ECDC4', 'd': '#45B7D1', 'f': '#FFA07A'}

# Plot 1: Quantum number relationships (top left)
ax1 = fig.add_subplot(gs[0, 0])
n_max = 4
data_points = []

for n in range(1, n_max + 1):
    for l in range(n):
        for ml in range(-l, l + 1):
            data_points.append((n, l, ml))

n_vals = [p[0] for p in data_points]
l_vals = [p[1] for p in data_points]
ml_vals = [p[2] for p in data_points]

# Create scatter plot
scatter = ax1.scatter(n_vals, l_vals, s=200, c=ml_vals, cmap='RdYlBu', 
                     edgecolors='black', linewidth=2, alpha=0.8)

ax1.set_xlabel('Principal quantum number (n)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Azimuthal quantum number (l)', fontsize=11, fontweight='bold')
ax1.set_title('Allowed Quantum Number Combinations', fontsize=13, fontweight='bold')
ax1.set_xticks(range(1, n_max + 1))
ax1.set_yticks(range(n_max))
ax1.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax1)
cbar.set_label('m_l value', fontsize=10)

# Plot 2: Orbital capacity (top middle)
ax2 = fig.add_subplot(gs[0, 1])
orbitals = ['s', 'p', 'd', 'f', 'g']
capacities = [2, 6, 10, 14, 18]

bars = ax2.bar(orbitals, capacities, color=[colors_orbital.get(o, '#95A5A6') for o in orbitals],
               edgecolor='black', linewidth=2, alpha=0.8)

ax2.set_xlabel('Orbital type', fontsize=11, fontweight='bold')
ax2.set_ylabel('Maximum electrons', fontsize=11, fontweight='bold')
ax2.set_title('Electron Capacity by Orbital Type', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, cap in zip(bars, capacities):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{cap}', ha='center', va='bottom', fontweight='bold', fontsize=11)

# Plot 3: Electron configuration progression (top right)
ax3 = fig.add_subplot(gs[0, 2])
elements_plot = list(range(1, 21))
total_electrons = elements_plot
s_electrons = []
p_electrons = []
d_electrons = []

for Z in elements_plot:
    config = generate_electron_configuration(Z)
    s_e = sum(count for orb, count in config.items() if 's' in orb)
    p_e = sum(count for orb, count in config.items() if 'p' in orb)
    d_e = sum(count for orb, count in config.items() if 'd' in orb)
    s_electrons.append(s_e)
    p_electrons.append(p_e)
    d_electrons.append(d_e)

ax3.fill_between(elements_plot, 0, s_electrons, label='s electrons', 
                color=colors_orbital['s'], alpha=0.7)
ax3.fill_between(elements_plot, s_electrons, 
                [s+p for s, p in zip(s_electrons, p_electrons)],
                label='p electrons', color=colors_orbital['p'], alpha=0.7)
ax3.fill_between(elements_plot, [s+p for s, p in zip(s_electrons, p_electrons)],
                [s+p+d for s, p, d in zip(s_electrons, p_electrons, d_electrons)],
                label='d electrons', color=colors_orbital['d'], alpha=0.7)

ax3.set_xlabel('Atomic Number (Z)', fontsize=11, fontweight='bold')
ax3.set_ylabel('Number of Electrons', fontsize=11, fontweight='bold')
ax3.set_title('Electron Configuration Build-up', fontsize=13, fontweight='bold')
ax3.legend(fontsize=10, loc='upper left')
ax3.grid(True, alpha=0.3)

# Plot 4: Energy level diagram (middle left)
ax4 = fig.add_subplot(gs[1, 0])
n_levels = 5
y_positions = [energy_level_hydrogen(n) for n in range(1, n_levels + 1)]

for i, (n, energy) in enumerate(zip(range(1, n_levels + 1), y_positions)):
    # Draw energy level line
    ax4.hlines(energy, n - 0.4, n + 0.4, colors='black', linewidth=3)
    
    # Label with quantum number and energy
    ax4.text(n + 0.5, energy, f'n={n}\n{energy:.2f} eV', 
            fontsize=9, va='center')
    
    # Draw orbitals
    orbitals_at_n = [orbital_name(n, l) for l in range(n)]
    for j, orb in enumerate(orbitals_at_n):
        x_offset = -0.3 + j * 0.2
        ax4.plot(n + x_offset, energy, 'o', markersize=10, 
                color=colors_orbital.get(orb[1], '#95A5A6'),
                markeredgecolor='black', markeredgewidth=1.5)

ax4.set_xlabel('Principal Quantum Number (n)', fontsize=11, fontweight='bold')
ax4.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax4.set_title('Hydrogen Atom Energy Levels', fontsize=13, fontweight='bold')
ax4.set_xticks(range(1, n_levels + 1))
ax4.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Ionization')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

# Plot 5: Radial probability distribution (middle middle)
ax5 = fig.add_subplot(gs[1, 1])
r = np.linspace(0, 30, 500)

orbitals_to_plot = [(1, 0), (2, 0), (2, 1), (3, 0)]
labels = ['1s', '2s', '2p', '3s']
colors_rad = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12']

for (n, l), label, color in zip(orbitals_to_plot, labels, colors_rad):
    R = radial_wavefunction(r, n, l)
    P = 4 * np.pi * r**2 * R**2  # Radial probability density
    ax5.plot(r, P, label=label, linewidth=2.5, color=color)

ax5.set_xlabel('Distance from nucleus (a₀)', fontsize=11, fontweight='bold')
ax5.set_ylabel('Radial Probability Density', fontsize=11, fontweight='bold')
ax5.set_title('Radial Probability Distributions', fontsize=13, fontweight='bold')
ax5.legend(fontsize=10)
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, 20)

# Plot 6: Hund's rule demonstration for Carbon (middle right)
ax6 = fig.add_subplot(gs[1, 2])
ax6.text(0.5, 0.95, "Carbon (Z=6): Hund's Rule", 
        ha='center', va='top', fontsize=13, fontweight='bold',
        transform=ax6.transAxes)

# Draw orbital boxes
orbital_sets = [
    ("1s", 1, [(0.5, 0.5), (0.5, -0.5)]),
    ("2s", 2, [(0.5, 0.5), (0.5, -0.5)]),
    ("2p", 3, [(0.5, 0), (0.5, 0)])  # Two electrons, parallel spins
]

y_base = 0.7
x_start = 0.25
box_width = 0.08
box_height = 0.15
arrow_size = 0.05

for orbital, num_boxes, electrons in orbital_sets:
    # Label
    ax6.text(0.1, y_base, orbital, fontsize=12, fontweight='bold',
            transform=ax6.transAxes, va='center')
    
    # Draw boxes
    for i in range(num_boxes):
        x = x_start + i * (box_width + 0.02)
        rect = Rectangle((x, y_base - box_height/2), box_width, box_height,
                            fill=False, edgecolor='black', linewidth=2,
                            transform=ax6.transAxes)
        ax6.add_patch(rect)
        ax6.add_patch(rect)
    
    # Draw electrons (arrows)
    for i, (spin_up, spin_down) in enumerate(electrons):
        if i >= num_boxes:
            break
        x = x_start + i * (box_width + 0.02) + box_width/2
        
        # Spin up
        if spin_up:
            ax6.arrow(x, y_base - 0.03, 0, arrow_size, 
                     head_width=0.02, head_length=0.015, fc='blue', ec='blue',
                     transform=ax6.transAxes)
        # Spin down
        if spin_down:
            ax6.arrow(x + 0.02, y_base + 0.03, 0, -arrow_size,
                     head_width=0.02, head_length=0.015, fc='red', ec='red',
                     transform=ax6.transAxes)
    
    y_base -= 0.25

ax6.text(0.5, 0.1, 
        "Hund's Rule: Electrons occupy\ndifferent orbitals with parallel\nspins before pairing",
        ha='center', va='center', fontsize=10, style='italic',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
        transform=ax6.transAxes)

ax6.axis('off')

# Plot 7: 3D orbital shapes - s orbital (bottom left)
ax7 = fig.add_subplot(gs[2, 0], projection='3d')

u = np.linspace(0, 2 * np.pi, 50)
v = np.linspace(0, np.pi, 50)
x_s = np.outer(np.cos(u), np.sin(v))
y_s = np.outer(np.sin(u), np.sin(v))
z_s = np.outer(np.ones(np.size(u)), np.cos(v))

ax7.plot_surface(x_s, y_s, z_s, color='#FF6B6B', alpha=0.8, edgecolor='none')
ax7.set_xlabel('x', fontsize=9)
ax7.set_ylabel('y', fontsize=9)
ax7.set_zlabel('z', fontsize=9)
ax7.set_title('s Orbital (l=0)\nSpherical', fontsize=12, fontweight='bold')
ax7.set_box_aspect([1,1,1])

# Plot 8: p orbital shape (bottom middle)
ax8 = fig.add_subplot(gs[2, 1], projection='3d')

u = np.linspace(0, 2 * np.pi, 50)
v = np.linspace(0, np.pi, 50)
# p_z orbital (dumbbell along z-axis)
x_p = 0.7 * np.outer(np.cos(u), np.sin(v)) * np.sin(v)
y_p = 0.7 * np.outer(np.sin(u), np.sin(v)) * np.sin(v)
z_p = 1.4 * np.outer(np.ones(np.size(u)), np.cos(v)**2) * np.sign(np.outer(np.ones(np.size(u)), np.cos(v)))

ax8.plot_surface(x_p, y_p, z_p, color='#4ECDC4', alpha=0.8, edgecolor='none')
ax8.plot_surface(x_p, y_p, -z_p, color='#4ECDC4', alpha=0.8, edgecolor='none')
ax8.set_xlabel('x', fontsize=9)
ax8.set_ylabel('y', fontsize=9)
ax8.set_zlabel('z', fontsize=9)
ax8.set_title('p Orbital (l=1)\nDumbbell', fontsize=12, fontweight='bold')
ax8.set_box_aspect([1,1,1])

plt.suptitle('Quantum Numbers and Electron States in Atoms', 
            fontsize=16, fontweight='bold', y=0.995)

plt.show(block=True)

print("\n" + "=" * 80)
print("KEY PRINCIPLES SUMMARY")
print("=" * 80)
print()
print("1. PAULI EXCLUSION PRINCIPLE:")
print("   • No two electrons can have identical quantum numbers")
print("   • Each orbital (n, l, m_l) holds maximum 2 electrons with opposite spins")
print()
print("2. AUFBAU PRINCIPLE:")
print("   • Fill lowest energy orbitals first")
print("   • Order: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d...")
print("   • Based on (n+l) rule: lower (n+l) fills first")
print()
print("3. HUND'S RULE:")
print("   • Maximize unpaired electrons in degenerate orbitals")
print("   • Electrons occupy different orbitals with parallel spins before pairing")
print("   • Minimizes electron-electron repulsion")
print()
print("4. QUANTUM NUMBER RULES:")
print("   • n ≥ 1 (principal quantum number)")
print("   • 0 ≤ l < n (azimuthal quantum number)")
print("   • -l ≤ m_l ≤ l (magnetic quantum number)")
print("   • m_s = ±½ (spin quantum number)")
print()
print("=" * 80)
print("KEY INSIGHTS:")
print("=" * 80)
print("• Each electron has a unique set of four quantum numbers (n, l, m_l, m_s)")
print("• Quantum numbers determine orbital energy, shape, and orientation")
print("• Pauli exclusion principle: no two electrons can have identical quantum numbers")
print("• Aufbau principle: electrons fill lowest energy orbitals first")
print("• Hund's rule: maximize unpaired spins in degenerate orbitals")
print("• s orbitals are spherical, p orbitals are dumbbell-shaped")
print("=" * 80)


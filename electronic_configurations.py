"""
Electronic Configurations: Aufbau Principle, Hund's Rule, and Pauli Exclusion
=============================================================================

This script demonstrates the fundamental principles governing electron
configurations in atoms, including:
- Aufbau principle (building-up principle)
- Hund's rule of maximum multiplicity
- Pauli exclusion principle

Introduction:
Electronic configuration describes the distribution of electrons among atomic
orbitals in an atom's ground state. Understanding these configurations is crucial
for predicting chemical properties, bonding behavior, and spectroscopic features.

Theory:

1. PAULI EXCLUSION PRINCIPLE (Wolfgang Pauli, 1925):
   No two electrons in an atom can have the same set of four quantum numbers
   (n, l, m_l, m_s). This fundamental principle arises from the antisymmetric
   nature of fermionic wavefunctions.
   
   Mathematical basis: For fermions, Ψ(r₁,r₂) = -Ψ(r₂,r₁)
   If two electrons had identical quantum numbers: Ψ = -Ψ → Ψ = 0 (impossible)
   
   Consequences:
   - Each orbital (defined by n, l, m_l) can hold maximum 2 electrons
   - These 2 electrons must have opposite spins (m_s = +½ and -½)
   - Maximum electrons in subshell: 2(2l+1)
     * s: 2 electrons (l=0, 1 orbital)
     * p: 6 electrons (l=1, 3 orbitals)
     * d: 10 electrons (l=2, 5 orbitals)
     * f: 14 electrons (l=3, 7 orbitals)

2. AUFBAU PRINCIPLE (from German "Aufbauprinzip" = building-up principle):
   Electrons occupy orbitals in order of increasing energy, filling lower
   energy orbitals before higher ones.
   
   Energy ordering rules:
   - Primary: Orbitals with lower (n + l) fill first (Madelung rule)
   - Secondary: If (n + l) is equal, lower n fills first
   
   Filling order: 1s → 2s → 2p → 3s → 3p → 4s → 3d → 4p → 5s → 4d → 5p → 6s → 4f → 5d → 6p → 7s...
   
   Why 4s before 3d? Although 3d is in a lower shell, (3+2)=5 > (4+0)=4
   The 4s orbital is more diffuse and has lower energy due to less electron-electron
   repulsion. However, once 3d starts filling, the actual energies become complex
   due to screening and exchange effects.
   
   Physical origin:
   - Nuclear-electron attraction (favors lower n)
   - Electron-electron repulsion (screening/shielding)
   - Penetration effect (s orbitals penetrate closer to nucleus)
   - Relativistic effects (important for heavy elements)

3. HUND'S RULE (Friedrich Hund, 1925):
   When filling degenerate orbitals (same energy), electrons occupy different
   orbitals with parallel spins before pairing.
   
   Full statement (three rules):
   
   Rule 1 (Maximum Multiplicity): For a given electron configuration, the term
   with maximum total spin S has the lowest energy. Electrons occupy orbitals
   singly with parallel spins before any pairing occurs.
   
   Physical basis:
   - Exchange interaction: Electrons with parallel spins stay farther apart
     spatially (Pauli principle for spatial wavefunction), reducing repulsion
   - Exchange energy: K ∝ ∫∫ ψᵢ*(r₁)ψⱼ*(r₂) 1/r₁₂ ψⱼ(r₁)ψᵢ(r₂) dr₁dr₂
   - Parallel spins → more exchange stabilization
   
   Rule 2 (Maximum Angular Momentum): For given multiplicity, the term with
   largest total orbital angular momentum L has lowest energy.
   
   Rule 3 (J-value): For given S and L:
   - If subshell less than half-filled: minimum J = |L - S| is lowest
   - If subshell more than half-filled: maximum J = L + S is lowest
   
   Examples:
   - Carbon (2p²): ↑ ↑ _ (triplet ³P) NOT ↑↓ _ (singlet ¹D or ¹S)
   - Nitrogen (2p³): ↑ ↑ ↑ (quartet ⁴S) - half-filled subshell extra stable
   - Oxygen (2p⁴): ↑↓ ↑ ↑ (triplet ³P) - start pairing from left

IMPORTANT CONSEQUENCES AND APPLICATIONS:

Electronic Configuration Principles Summary:
   • Pauli Exclusion: Maximum 2 electrons per orbital (opposite spins: ↑↓)
   • Aufbau: Fill orbitals in order of increasing (n+l); if equal, lower n first
   • Hund's Rule: Maximize unpaired electrons in degenerate orbitals (parallel spins)

Special Stability Cases:
   • Half-filled subshells: N (2p³), P (3p³), Cr (3d⁵), Mn (3d⁵)
     - All spins parallel, maximum exchange stabilization
   • Fully-filled subshells: He (1s²), Ne ([He]2s²2p⁶), Cu (3d¹⁰), Zn (3d¹⁰)
     - Closed shell, extra stability
   • Noble gases: He, Ne, Ar, Kr, Xe, Rn - filled shells, chemically inert

Exceptions to Aufbau:
   • Chromium (Z=24): [Ar] 4s¹ 3d⁵ instead of [Ar] 4s² 3d⁴
     Reason: Half-filled d⁵ provides extra stability
   • Copper (Z=29): [Ar] 4s¹ 3d¹⁰ instead of [Ar] 4s² 3d⁹
     Reason: Fully-filled d¹⁰ provides extra stability

Magnetic Properties:
   • Paramagnetic: Unpaired electrons → attracted to magnetic field
     Examples: O₂ (2 unpaired), Fe (4 unpaired), Cu²⁺ (1 unpaired)
   • Diamagnetic: All electrons paired → slightly repelled by magnetic field
     Examples: N₂, noble gases, Zn

Spin Multiplicity:
   • Multiplicity = 2S + 1, where S = total spin quantum number
   • Singlet (S=0): All electrons paired (He, Be, Ne)
   • Doublet (S=½): 1 unpaired electron (H, Li, Na, Cu)
   • Triplet (S=1): 2 unpaired electrons (C, O, S)
   • Quartet (S=3/2): 3 unpaired electrons (N, P)
   • Higher multiplicities in transition metals with many d electrons

Chemical Reactivity:
   • Valence electrons determine bonding behavior
   • Elements with similar valence configurations have similar chemistry
   • Transition metals: d electrons participate in bonding, variable oxidation states
   • Lanthanides/Actinides: f electrons shield poorly, similar chemistry within series

Author: ewdhp
Date: November 11, 2025
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Interactive backend
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

# ============================================================================
# CONSTANTS AND DATA
# ============================================================================

# Orbital energy ordering (simplified, for building configurations)
ORBITAL_ORDER = [
    '1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s',
    '4f', '5d', '6p', '7s', '5f', '6d', '7p'
]

# Maximum electrons per orbital type
MAX_ELECTRONS = {
    's': 2, 'p': 6, 'd': 10, 'f': 14
}

# Element data (Z: [symbol, name])
ELEMENTS = {
    1: ['H', 'Hydrogen'], 2: ['He', 'Helium'], 3: ['Li', 'Lithium'],
    4: ['Be', 'Beryllium'], 5: ['B', 'Boron'], 6: ['C', 'Carbon'],
    7: ['N', 'Nitrogen'], 8: ['O', 'Oxygen'], 9: ['F', 'Fluorine'],
    10: ['Ne', 'Neon'], 11: ['Na', 'Sodium'], 12: ['Mg', 'Magnesium'],
    13: ['Al', 'Aluminum'], 14: ['Si', 'Silicon'], 15: ['P', 'Phosphorus'],
    16: ['S', 'Sulfur'], 17: ['Cl', 'Chlorine'], 18: ['Ar', 'Argon'],
    19: ['K', 'Potassium'], 20: ['Ca', 'Calcium'], 21: ['Sc', 'Scandium'],
    22: ['Ti', 'Titanium'], 23: ['V', 'Vanadium'], 24: ['Cr', 'Chromium'],
    25: ['Mn', 'Manganese'], 26: ['Fe', 'Iron'], 27: ['Co', 'Cobalt'],
    28: ['Ni', 'Nickel'], 29: ['Cu', 'Copper'], 30: ['Zn', 'Zinc'],
    31: ['Ga', 'Gallium'], 32: ['Ge', 'Germanium'], 33: ['As', 'Arsenic'],
    34: ['Se', 'Selenium'], 35: ['Br', 'Bromine'], 36: ['Kr', 'Krypton']
}

# Noble gas cores
NOBLE_GAS_CORES = {
    'He': 2, 'Ne': 10, 'Ar': 18, 'Kr': 36, 'Xe': 54, 'Rn': 86
}

# ============================================================================
# FUNCTIONS
# ============================================================================

def get_orbital_type(orbital_str):
    """Extract orbital type (s, p, d, f) from orbital string like '3d'."""
    return orbital_str[-1]

def get_n_quantum_number(orbital_str):
    """Extract principal quantum number from orbital string."""
    return int(orbital_str[:-1])

def build_electron_configuration(Z):
    """
    Build ground state electron configuration for element with atomic number Z.
    
    Parameters:
    -----------
    Z : int
        Atomic number (number of electrons)
    
    Returns:
    --------
    config : dict
        Dictionary mapping orbital names to electron counts
    """
    config = {}
    electrons_remaining = Z
    
    for orbital in ORBITAL_ORDER:
        if electrons_remaining <= 0:
            break
        
        orbital_type = get_orbital_type(orbital)
        max_e = MAX_ELECTRONS[orbital_type]
        
        # Fill orbital with min(electrons_remaining, max_capacity)
        electrons_in_orbital = min(electrons_remaining, max_e)
        config[orbital] = electrons_in_orbital
        electrons_remaining -= electrons_in_orbital
    
    return config

def config_to_string(config, use_superscript=True):
    """
    Convert configuration dictionary to standard notation.
    
    Parameters:
    -----------
    config : dict
        Electron configuration
    use_superscript : bool
        If True, use superscript numbers (e.g., 2s²)
    
    Returns:
    --------
    str : Configuration string
    """
    if use_superscript:
        superscripts = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        parts = [f"{orb}{count}".translate(superscripts) 
                for orb, count in config.items() if count > 0]
    else:
        parts = [f"{orb}^{count}" for orb, count in config.items() if count > 0]
    
    return " ".join(parts)

def config_to_noble_gas_notation(config, Z):
    """
    Convert to noble gas core notation (e.g., [Ne] 3s² 3p⁴ for S).
    
    Parameters:
    -----------
    config : dict
        Electron configuration
    Z : int
        Atomic number
    
    Returns:
    --------
    str : Noble gas notation
    """
    # Find the largest noble gas with Z_ng < Z
    noble_gas = None
    noble_z = 0
    
    for gas, z in sorted(NOBLE_GAS_CORES.items(), key=lambda x: x[1], reverse=True):
        if z < Z:
            noble_gas = gas
            noble_z = z
            break
    
    if noble_gas is None:
        return config_to_string(config)
    
    # Build configuration for remaining electrons
    remaining_config = build_electron_configuration(Z - noble_z)
    
    return f"[{noble_gas}] {config_to_string(remaining_config)}"

def get_valence_subshell(config):
    """
    Identify the outermost partially filled subshell.
    
    Returns:
    --------
    tuple : (orbital_name, n_electrons, max_electrons)
    """
    # Find last orbital with electrons
    last_orbital = None
    for orbital in reversed(ORBITAL_ORDER):
        if orbital in config and config[orbital] > 0:
            orbital_type = get_orbital_type(orbital)
            max_e = MAX_ELECTRONS[orbital_type]
            
            # Return first partially filled or last filled
            if config[orbital] < max_e or last_orbital is None:
                return (orbital, config[orbital], max_e)
            last_orbital = orbital
    
    return (last_orbital, config[last_orbital], MAX_ELECTRONS[get_orbital_type(last_orbital)])

def get_orbital_box_diagram(orbital, n_electrons):
    """
    Generate orbital box diagram following Hund's rule.
    
    Parameters:
    -----------
    orbital : str
        Orbital name (e.g., '2p', '3d')
    n_electrons : int
        Number of electrons in this orbital
    
    Returns:
    --------
    list : List of tuples (n_up, n_down) for each orbital
    """
    orbital_type = get_orbital_type(orbital)
    
    # Number of spatial orbitals
    n_orbitals = {
        's': 1,  # m_l = 0
        'p': 3,  # m_l = -1, 0, +1
        'd': 5,  # m_l = -2, -1, 0, +1, +2
        'f': 7   # m_l = -3, -2, -1, 0, +1, +2, +3
    }[orbital_type]
    
    # Initialize empty orbitals
    boxes = [[0, 0] for _ in range(n_orbitals)]  # [spin_up, spin_down]
    
    electrons_remaining = n_electrons
    
    # Hund's Rule 1: Fill all orbitals singly with parallel spins first
    for i in range(n_orbitals):
        if electrons_remaining > 0:
            boxes[i][0] = 1  # Spin up
            electrons_remaining -= 1
        else:
            break
    
    # Now pair electrons (spin down)
    for i in range(n_orbitals):
        if electrons_remaining > 0:
            boxes[i][1] = 1  # Spin down
            electrons_remaining -= 1
        else:
            break
    
    return boxes

def calculate_total_spin(config):
    """
    Calculate total spin quantum number S for ground state.
    
    Returns:
    --------
    float : Total spin S (in units of ℏ)
    """
    unpaired = 0
    
    for orbital, n_electrons in config.items():
        orbital_type = get_orbital_type(orbital)
        max_e = MAX_ELECTRONS[orbital_type]
        n_orbitals = max_e // 2
        
        # Count unpaired electrons
        if n_electrons <= n_orbitals:
            # All electrons unpaired (Hund's rule)
            unpaired += n_electrons
        else:
            # Some paired, some unpaired
            unpaired += max_e - n_electrons
    
    return unpaired / 2.0

def calculate_multiplicity(S):
    """
    Calculate spin multiplicity: 2S + 1.
    
    Parameters:
    -----------
    S : float
        Total spin quantum number
    
    Returns:
    --------
    int : Multiplicity
    """
    return int(2 * S + 1)

def is_noble_gas(Z):
    """Check if element is a noble gas."""
    return Z in NOBLE_GAS_CORES.values()

def is_half_filled_or_filled(config, orbital):
    """
    Check if orbital is half-filled or fully filled (extra stability).
    """
    if orbital not in config:
        return False
    
    orbital_type = get_orbital_type(orbital)
    max_e = MAX_ELECTRONS[orbital_type]
    n_e = config[orbital]
    
    return n_e == max_e // 2 or n_e == max_e

# ============================================================================
# MAIN DEMONSTRATION
# ============================================================================

print("=" * 80)
print("ELECTRONIC CONFIGURATIONS: AUFBAU, HUND'S RULE, AND PAULI EXCLUSION")
print("=" * 80)
print("\nDemonstrating the fundamental principles governing electron configurations")
print("in ground state atoms.\n")

# ============================================================================
# PART 1: ORBITAL FILLING ORDER (AUFBAU PRINCIPLE)
# ============================================================================

print("1. AUFBAU PRINCIPLE - ORBITAL FILLING ORDER")
print("-" * 80)
print("\nThe (n + l) rule determines orbital filling order:")
print()
print(f"{'Orbital':<10} {'n':<5} {'l':<5} {'n+l':<8} {'Order':<8} {'Max e⁻':<10}")
print("-" * 80)

order_data = []
for i, orbital in enumerate(ORBITAL_ORDER[:15], 1):  # First 15 orbitals
    n = get_n_quantum_number(orbital)
    l_values = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    l = l_values[get_orbital_type(orbital)]
    orbital_type = get_orbital_type(orbital)
    max_e = MAX_ELECTRONS[orbital_type]
    
    print(f"{orbital:<10} {n:<5} {l:<5} {n+l:<8} {i:<8} {max_e:<10}")
    order_data.append((orbital, n, l, n+l, i, max_e))

print("\nNote: Orbitals are filled in order of increasing (n+l).")
print("      When (n+l) is equal, lower n fills first.")
print()

# ============================================================================
# PART 2: ELECTRON CONFIGURATIONS FOR ELEMENTS
# ============================================================================

print("2. ELECTRON CONFIGURATIONS (Z = 1 to 36)")
print("-" * 80)
print()

# Build configurations for elements
configs_full = {}
configs_ng = {}

print(f"{'Z':<4} {'Symbol':<6} {'Name':<12} {'Configuration':<35} {'Noble Gas Notation':<30}")
print("-" * 80)

demo_elements = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 18, 19, 20, 24, 25, 26, 29, 30, 36]

for Z in demo_elements:
    if Z in ELEMENTS:
        symbol, name = ELEMENTS[Z]
        config = build_electron_configuration(Z)
        config_str = config_to_string(config)
        ng_notation = config_to_noble_gas_notation(config, Z)
        
        configs_full[Z] = config
        configs_ng[Z] = ng_notation
        
        print(f"{Z:<4} {symbol:<6} {name:<12} {config_str:<35} {ng_notation:<30}")

print()

# ============================================================================
# PART 3: SPECIAL CASES AND EXCEPTIONS
# ============================================================================

print("3. EXCEPTIONS TO AUFBAU PRINCIPLE")
print("-" * 80)
print("\nSome elements have anomalous configurations due to extra stability of")
print("half-filled and fully-filled subshells.\n")

exceptions = [
    (24, 'Cr', '[Ar] 4s¹ 3d⁵', '[Ar] 4s² 3d⁴', 'Half-filled d⁵ more stable'),
    (29, 'Cu', '[Ar] 4s¹ 3d¹⁰', '[Ar] 4s² 3d⁹', 'Filled d¹⁰ more stable'),
    (42, 'Mo', '[Kr] 5s¹ 4d⁵', '[Kr] 5s² 4d⁴', 'Half-filled d⁵ more stable'),
    (47, 'Ag', '[Kr] 5s¹ 4d¹⁰', '[Kr] 5s² 4d⁹', 'Filled d¹⁰ more stable'),
]

print(f"{'Z':<4} {'Element':<8} {'Actual Config':<20} {'Expected Config':<20} {'Reason':<30}")
print("-" * 80)
for Z, elem, actual, expected, reason in exceptions:
    print(f"{Z:<4} {elem:<8} {actual:<20} {expected:<20} {reason:<30}")

print()

# ============================================================================
# PART 4: HUND'S RULE EXAMPLES
# ============================================================================

print("4. HUND'S RULE - ORBITAL BOX DIAGRAMS")
print("-" * 80)
print("\nElectrons occupy degenerate orbitals singly with parallel spins")
print("before pairing (maximizes total spin).\n")

hunds_examples = [
    (5, 'B', '2p¹'),
    (6, 'C', '2p²'),
    (7, 'N', '2p³'),
    (8, 'O', '2p⁴'),
    (9, 'F', '2p⁵'),
    (15, 'P', '3p³'),
    (16, 'S', '3p⁴'),
    (26, 'Fe', '3d⁶'),
]

for Z, symbol, valence in hunds_examples:
    config = build_electron_configuration(Z)
    orbital, n_e, max_e = get_valence_subshell(config)
    boxes = get_orbital_box_diagram(orbital, n_e)
    S = calculate_total_spin(config)
    multiplicity = calculate_multiplicity(S)
    
    # Create box diagram string
    box_str = ""
    for up, down in boxes:
        arrows = ""
        if up and down:
            arrows = "↑↓"
        elif up:
            arrows = "↑ "
        elif down:
            arrows = " ↓"
        else:
            arrows = "  "
        box_str += f"[{arrows}] "
    
    print(f"{symbol} (Z={Z}): {orbital:<4} {box_str}  S={S:.1f}, Multiplicity={multiplicity}")

print()

# ============================================================================
# PART 5: PAULI EXCLUSION DEMONSTRATION
# ============================================================================

print("5. PAULI EXCLUSION PRINCIPLE")
print("-" * 80)
print("\nMaximum 2 electrons per orbital (with opposite spins).\n")

print("Maximum electrons in each subshell:")
print(f"{'Subshell':<12} {'l value':<10} {'# orbitals':<15} {'Max e⁻ (2×orbitals)':<20}")
print("-" * 80)

for orbital_type in ['s', 'p', 'd', 'f']:
    l = {'s': 0, 'p': 1, 'd': 2, 'f': 3}[orbital_type]
    n_orbitals = 2 * l + 1
    max_e = 2 * n_orbitals
    print(f"{orbital_type:<12} {l:<10} {n_orbitals:<15} {max_e:<20}")

print("\nEach orbital can hold max 2 electrons with opposite spins:")
print("  Spin up: m_s = +½")
print("  Spin down: m_s = -½")
print()

# ============================================================================
# VISUALIZATIONS
# ============================================================================

print("\nGenerating visualizations...")
print("=" * 80)

fig = plt.figure(figsize=(20, 12))
gs = GridSpec(3, 4, figure=fig, hspace=0.4, wspace=0.4)

# Color scheme
color_up = '#3498DB'    # Blue for spin up
color_down = '#E74C3C'  # Red for spin down
color_orbital_s = '#FF6B6B'
color_orbital_p = '#4ECDC4'
color_orbital_d = '#45B7D1'
color_orbital_f = '#FFA07A'

# ============================================================================
# PLOT 1: AUFBAU PRINCIPLE - ORBITAL ENERGY DIAGRAM
# ============================================================================

ax1 = fig.add_subplot(gs[0:2, 0])

# Energy levels (approximate, for visualization)
energy_levels = {
    '1s': 1, '2s': 2, '2p': 2.5, '3s': 3, '3p': 3.5,
    '4s': 4, '3d': 4.2, '4p': 4.5, '5s': 5, '4d': 5.2,
    '5p': 5.5, '6s': 6, '4f': 6.1, '5d': 6.3, '6p': 6.6
}

y_pos = 0
orbital_positions = {}

for orbital in ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d']:
    if orbital in energy_levels:
        E = energy_levels[orbital]
        orbital_type = get_orbital_type(orbital)
        
        # Determine color
        colors = {'s': color_orbital_s, 'p': color_orbital_p, 
                 'd': color_orbital_d, 'f': color_orbital_f}
        color = colors[orbital_type]
        
        # Draw energy level line
        ax1.plot([0, 1], [E, E], linewidth=3, color=color, alpha=0.8)
        
        # Label
        ax1.text(1.1, E, orbital, fontsize=11, fontweight='bold', va='center')
        
        orbital_positions[orbital] = E

# Draw filling arrows
arrow_x = -0.3
for i, orbital in enumerate(['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p']):
    if i > 0 and orbital in energy_levels:
        prev_orbital = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p'][i-1]
        if prev_orbital in orbital_positions:
            y1 = orbital_positions[prev_orbital]
            y2 = orbital_positions[orbital]
            ax1.annotate('', xy=(arrow_x, y2), xytext=(arrow_x, y1),
                        arrowprops=dict(arrowstyle='->', lw=1.5, color='gray', alpha=0.6))

ax1.set_xlim(-0.5, 1.8)
ax1.set_ylim(0.5, 7)
ax1.set_ylabel('Energy (arbitrary units)', fontsize=12, fontweight='bold')
ax1.set_title('Aufbau Principle:\nOrbital Filling Order', fontsize=13, fontweight='bold')
ax1.set_xticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

# ============================================================================
# PLOT 2: (n+l) RULE VISUALIZATION
# ============================================================================

ax2 = fig.add_subplot(gs[0:2, 1])

# Create grid showing (n, l) combinations
n_max = 6
l_max = 3

for n in range(1, n_max + 1):
    for l in range(min(n, l_max + 1)):
        nl_sum = n + l
        
        # Position
        x = n
        y = l
        
        # Color by (n+l) value
        orbital_type = ['s', 'p', 'd', 'f'][l]
        orbital_name = f"{n}{orbital_type}"
        
        # Draw box
        colors_nl = {2: '#FFCDD2', 3: '#F8BBD0', 4: '#E1BEE7', 5: '#D1C4E9',
                    6: '#C5CAE9', 7: '#BBDEFB', 8: '#B3E5FC'}
        color = colors_nl.get(nl_sum, '#E0E0E0')
        
        rect = Rectangle((x - 0.4, y - 0.4), 0.8, 0.8, 
                         facecolor=color, edgecolor='black', linewidth=2)
        ax2.add_patch(rect)
        
        # Label
        ax2.text(x, y, orbital_name, ha='center', va='center', 
                fontsize=12, fontweight='bold')
        
        # Show (n+l) value
        ax2.text(x, y - 0.28, f"n+l={nl_sum}", ha='center', va='top',
                fontsize=8, style='italic')

ax2.set_xlim(0.5, n_max + 0.5)
ax2.set_ylim(-0.5, l_max + 0.5)
ax2.set_xlabel('Principal quantum number (n)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Angular momentum (l)', fontsize=11, fontweight='bold')
ax2.set_yticks(range(l_max + 1))
ax2.set_yticklabels(['0 (s)', '1 (p)', '2 (d)', '3 (f)'])
ax2.set_xticks(range(1, n_max + 1))
ax2.set_title('(n + l) Rule:\nOrbitals with Lower (n+l) Fill First', 
             fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)

# ============================================================================
# PLOT 3: HUND'S RULE - CARBON EXAMPLE
# ============================================================================

ax3 = fig.add_subplot(gs[0, 2])
ax3.text(0.5, 0.95, "Hund's Rule: Carbon (2p²)", 
        ha='center', va='top', fontsize=13, fontweight='bold',
        transform=ax3.transAxes)

# Show two possible configurations
y_correct = 0.65
y_wrong = 0.25

# Correct (Hund's rule)
ax3.text(0.05, y_correct + 0.05, "✓ CORRECT (³P triplet)", 
        transform=ax3.transAxes, fontsize=11, fontweight='bold', color='green')

boxes_correct = [[1, 0], [1, 0], [0, 0]]  # ↑ ↑ _
x_start = 0.2
for i, (up, down) in enumerate(boxes_correct):
    x = x_start + i * 0.22
    
    # Draw box
    rect = Rectangle((x, y_correct - 0.05), 0.15, 0.1,
                     facecolor='white', edgecolor='black', linewidth=2,
                     transform=ax3.transAxes)
    ax3.add_patch(rect)
    
    # Draw arrows
    if up:
        ax3.arrow(x + 0.075, y_correct - 0.02, 0, 0.04,
                 head_width=0.025, head_length=0.015, fc=color_up, ec=color_up,
                 transform=ax3.transAxes)
    if down:
        ax3.arrow(x + 0.075, y_correct + 0.02, 0, -0.04,
                 head_width=0.025, head_length=0.015, fc=color_down, ec=color_down,
                 transform=ax3.transAxes)

ax3.text(0.87, y_correct, "S = 1, Multiplicity = 3",
        transform=ax3.transAxes, fontsize=10, va='center')

# Incorrect (paired)
ax3.text(0.05, y_wrong + 0.05, "✗ WRONG (¹D singlet)", 
        transform=ax3.transAxes, fontsize=11, fontweight='bold', color='red')

boxes_wrong = [[1, 1], [0, 0], [0, 0]]  # ↑↓ _ _
for i, (up, down) in enumerate(boxes_wrong):
    x = x_start + i * 0.22
    
    # Draw box
    rect = Rectangle((x, y_wrong - 0.05), 0.15, 0.1,
                     facecolor='white', edgecolor='black', linewidth=2,
                     transform=ax3.transAxes)
    ax3.add_patch(rect)
    
    # Draw arrows
    if up:
        ax3.arrow(x + 0.075, y_wrong - 0.02, 0, 0.04,
                 head_width=0.025, head_length=0.015, fc=color_up, ec=color_up,
                 transform=ax3.transAxes)
    if down:
        ax3.arrow(x + 0.075, y_wrong + 0.02, 0, -0.04,
                 head_width=0.025, head_length=0.015, fc=color_down, ec=color_down,
                 transform=ax3.transAxes)

ax3.text(0.87, y_wrong, "S = 0, Multiplicity = 1",
        transform=ax3.transAxes, fontsize=10, va='center')

ax3.axis('off')

# ============================================================================
# PLOT 4: NITROGEN - HALF-FILLED STABILITY
# ============================================================================

ax4 = fig.add_subplot(gs[0, 3])
ax4.text(0.5, 0.95, "Nitrogen (2p³)\nHalf-Filled Stability", 
        ha='center', va='top', fontsize=13, fontweight='bold',
        transform=ax4.transAxes)

y_n = 0.5
boxes_n = [[1, 0], [1, 0], [1, 0]]  # ↑ ↑ ↑ (all parallel)
x_start = 0.15

for i, (up, down) in enumerate(boxes_n):
    x = x_start + i * 0.22
    
    # Draw box with special color for half-filled
    rect = Rectangle((x, y_n - 0.08), 0.15, 0.16,
                     facecolor='#FFFFCC', edgecolor='green', linewidth=3,
                     transform=ax4.transAxes)
    ax4.add_patch(rect)
    
    # Draw arrows
    if up:
        ax3.arrow(x + 0.075, y_n - 0.03, 0, 0.08,
                 head_width=0.025, head_length=0.02, fc=color_up, ec=color_up,
                 transform=ax4.transAxes)

ax4.text(0.5, y_n - 0.2, "All spins parallel\nS = 3/2, Multiplicity = 4 (⁴S)\nExtra stable!",
        ha='center', va='top', fontsize=10, style='italic',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5),
        transform=ax4.transAxes)

ax4.axis('off')

# ============================================================================
# PLOT 5: OXYGEN - START PAIRING
# ============================================================================

ax5 = fig.add_subplot(gs[1, 2])
ax5.text(0.5, 0.95, "Oxygen (2p⁴)\nStart Pairing", 
        ha='center', va='top', fontsize=13, fontweight='bold',
        transform=ax5.transAxes)

y_o = 0.5
boxes_o = [[1, 1], [1, 0], [1, 0]]  # ↑↓ ↑ ↑
x_start = 0.15

for i, (up, down) in enumerate(boxes_o):
    x = x_start + i * 0.22
    
    # Draw box
    rect = Rectangle((x, y_o - 0.08), 0.15, 0.16,
                     facecolor='white', edgecolor='black', linewidth=2,
                     transform=ax5.transAxes)
    ax5.add_patch(rect)
    
    # Draw arrows
    if up:
        ax5.arrow(x + 0.075, y_o - 0.03, 0, 0.08,
                 head_width=0.025, head_length=0.02, fc=color_up, ec=color_up,
                 transform=ax5.transAxes)
    if down:
        ax5.arrow(x + 0.075, y_o + 0.05, 0, -0.08,
                 head_width=0.025, head_length=0.02, fc=color_down, ec=color_down,
                 transform=ax5.transAxes)

ax5.text(0.5, y_o - 0.2, "4th electron pairs\nS = 1, Multiplicity = 3 (³P)",
        ha='center', va='top', fontsize=10, style='italic',
        transform=ax5.transAxes)

ax5.axis('off')

# ============================================================================
# PLOT 6: PAULI EXCLUSION - ORBITAL CAPACITY
# ============================================================================

ax6 = fig.add_subplot(gs[1, 3])
ax6.text(0.5, 0.95, "Pauli Exclusion Principle", 
        ha='center', va='top', fontsize=13, fontweight='bold',
        transform=ax6.transAxes)

# Show single orbital with max 2 electrons
y_pauli = 0.6
x_center = 0.5

# Draw orbital box
rect = Rectangle((x_center - 0.1, y_pauli - 0.08), 0.2, 0.16,
                 facecolor='white', edgecolor='purple', linewidth=3,
                 transform=ax6.transAxes)
ax6.add_patch(rect)

# Two electrons with opposite spins
ax6.arrow(x_center, y_pauli - 0.03, 0, 0.08,
         head_width=0.03, head_length=0.02, fc=color_up, ec=color_up,
         transform=ax6.transAxes)
ax6.arrow(x_center, y_pauli + 0.05, 0, -0.08,
         head_width=0.03, head_length=0.02, fc=color_down, ec=color_down,
         transform=ax6.transAxes)

ax6.text(0.5, y_pauli - 0.2, 
        "One orbital holds\nMAX 2 electrons\nwith opposite spins",
        ha='center', va='top', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.7),
        transform=ax6.transAxes)

ax6.text(0.5, 0.1, "m_s = +½ (↑) and m_s = -½ (↓)",
        ha='center', va='center', fontsize=9, style='italic',
        transform=ax6.transAxes)

ax6.axis('off')

# ============================================================================
# PLOT 7: CONFIGURATION BUILD-UP (H to Ne)
# ============================================================================

ax7 = fig.add_subplot(gs[2, :2])

elements_demo = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y_spacing = 1.0

ax7.text(0.5, 10.8, 'Electron Configuration Build-up (H → Ne)', 
        ha='center', fontsize=13, fontweight='bold')

for i, Z in enumerate(elements_demo):
    symbol, name = ELEMENTS[Z]
    config = build_electron_configuration(Z)
    config_str = config_to_string(config)
    
    y = 10 - i * y_spacing
    
    # Element label
    ax7.text(-1.5, y, f"{symbol} (Z={Z})", fontsize=10, fontweight='bold', va='center')
    
    # Draw orbital boxes
    x_pos = 0
    
    # 1s orbital
    if '1s' in config:
        n_e = config['1s']
        boxes_1s = get_orbital_box_diagram('1s', n_e)
        
        for j, (up, down) in enumerate(boxes_1s):
            rect = Rectangle((x_pos, y - 0.3), 0.6, 0.6,
                           facecolor=color_orbital_s, edgecolor='black', 
                           linewidth=1.5, alpha=0.3)
            ax7.add_patch(rect)
            
            if up:
                ax7.arrow(x_pos + 0.3, y - 0.1, 0, 0.3,
                         head_width=0.15, head_length=0.1, fc=color_up, ec=color_up)
            if down:
                ax7.arrow(x_pos + 0.3, y + 0.1, 0, -0.3,
                         head_width=0.15, head_length=0.1, fc=color_down, ec=color_down)
            
            if i == 0:  # Label for first element
                ax7.text(x_pos + 0.3, y - 0.5, '1s', ha='center', fontsize=8)
        
        x_pos += 1.0
    
    # 2s orbital
    if '2s' in config:
        n_e = config['2s']
        boxes_2s = get_orbital_box_diagram('2s', n_e)
        
        for j, (up, down) in enumerate(boxes_2s):
            rect = Rectangle((x_pos, y - 0.3), 0.6, 0.6,
                           facecolor=color_orbital_s, edgecolor='black', 
                           linewidth=1.5, alpha=0.3)
            ax7.add_patch(rect)
            
            if up:
                ax7.arrow(x_pos + 0.3, y - 0.1, 0, 0.3,
                         head_width=0.15, head_length=0.1, fc=color_up, ec=color_up)
            if down:
                ax7.arrow(x_pos + 0.3, y + 0.1, 0, -0.3,
                         head_width=0.15, head_length=0.1, fc=color_down, ec=color_down)
            
            if Z == 3:  # Label for Li
                ax7.text(x_pos + 0.3, y - 0.5, '2s', ha='center', fontsize=8)
        
        x_pos += 1.0
    
    # 2p orbitals
    if '2p' in config:
        n_e = config['2p']
        boxes_2p = get_orbital_box_diagram('2p', n_e)
        
        for j, (up, down) in enumerate(boxes_2p):
            rect = Rectangle((x_pos + j * 0.7, y - 0.3), 0.6, 0.6,
                           facecolor=color_orbital_p, edgecolor='black', 
                           linewidth=1.5, alpha=0.3)
            ax7.add_patch(rect)
            
            if up:
                ax7.arrow(x_pos + j * 0.7 + 0.3, y - 0.1, 0, 0.3,
                         head_width=0.15, head_length=0.1, fc=color_up, ec=color_up)
            if down:
                ax7.arrow(x_pos + j * 0.7 + 0.3, y + 0.1, 0, -0.3,
                         head_width=0.15, head_length=0.1, fc=color_down, ec=color_down)
            
            if Z == 5:  # Label for B
                ax7.text(x_pos + j * 0.7 + 0.3, y - 0.5, '2p', ha='center', fontsize=8)
    
    # Configuration string
    ax7.text(7.5, y, config_str, fontsize=10, va='center', family='monospace')

ax7.set_xlim(-2, 10)
ax7.set_ylim(-0.5, 11)
ax7.axis('off')

# ============================================================================
# PLOT 8: MULTIPLICITY vs ATOMIC NUMBER
# ============================================================================

ax8 = fig.add_subplot(gs[2, 2:])

Z_range = range(1, 37)
multiplicities = []

for Z in Z_range:
    config = build_electron_configuration(Z)
    S = calculate_total_spin(config)
    mult = calculate_multiplicity(S)
    multiplicities.append(mult)

ax8.plot(Z_range, multiplicities, 'o-', linewidth=2, markersize=6, 
        color='#9B59B6', markerfacecolor='#E8DAEF', markeredgewidth=2)

# Highlight noble gases
for Z in [2, 10, 18, 36]:
    if Z in Z_range:
        ax8.axvline(Z, color='gold', linestyle='--', linewidth=2, alpha=0.5)
        symbol = ELEMENTS[Z][0]
        ax8.text(Z, max(multiplicities) * 0.95, symbol, 
                ha='center', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='gold', alpha=0.3))

ax8.set_xlabel('Atomic Number (Z)', fontsize=11, fontweight='bold')
ax8.set_ylabel('Spin Multiplicity (2S+1)', fontsize=11, fontweight='bold')
ax8.set_title('Spin Multiplicity vs Atomic Number\n(Singlets, Doublets, Triplets, Quartets)', 
             fontsize=12, fontweight='bold')
ax8.grid(True, alpha=0.3)
ax8.set_xlim(0, 37)
ax8.set_ylim(0.5, max(multiplicities) + 0.5)

# Add multiplicity labels
mult_names = {1: 'Singlet', 2: 'Doublet', 3: 'Triplet', 4: 'Quartet', 5: 'Quintet', 6: 'Sextet', 7: 'Septet'}
for mult, name in mult_names.items():
    if mult <= max(multiplicities):
        ax8.text(36.5, mult, name, fontsize=9, va='center', style='italic')

plt.suptitle('Electronic Configurations: Aufbau Principle, Hund\'s Rule, and Pauli Exclusion', 
            fontsize=16, fontweight='bold', y=0.995)

plt.show(block=True)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("KEY INSIGHTS:")
print("=" * 80)
print("• Pauli Exclusion: Max 2 electrons per orbital with opposite spins")
print("• Aufbau Principle: Fill orbitals in order of increasing (n+l), then n")
print("• Hund's Rule: Maximize unpaired electrons in degenerate orbitals")
print("• Half-filled and fully-filled subshells have extra stability")
print("• Multiplicity = 2S+1, where S is total spin quantum number")
print("• Noble gases have closed shells (multiplicity = 1, diamagnetic)")
print("• Transition metals often have unpaired d electrons (paramagnetic)")
print("=" * 80)
print("\nVisualization complete!")
print("=" * 80)

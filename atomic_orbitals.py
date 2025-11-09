"""
Atomic Orbitals: Shapes and Energies (s, p, d, f)
=================================================

This script demonstrates the shapes and energies of atomic orbitals in hydrogen-like atoms.
We explore the complete wavefunction Ψ(r,θ,φ) = R(r)Y(θ,φ), including:
- Radial component R_nl(r)
- Angular component Y_lm(θ,φ) (spherical harmonics)
- 3D probability density distributions
- Nodal surfaces and orbital energies

Author: Educational Demonstration
License: MIT
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Use interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.special import sph_harm, genlaguerre, factorial
from matplotlib.gridspec import GridSpec

# Physical constants
BOHR_RADIUS = 5.29177210903e-11  # meters
RYDBERG_ENERGY = 13.605693122994  # eV
HBAR = 1.054571817e-34  # J·s
ELECTRON_MASS = 9.1093837015e-31  # kg

def radial_wavefunction(n, l, r, Z=1):
    """
    Calculate the radial wavefunction R_nl(r) for hydrogen-like atoms.
    
    Parameters:
    -----------
    n : int
        Principal quantum number (n ≥ 1)
    l : int
        Angular momentum quantum number (0 ≤ l < n)
    r : ndarray
        Radial distance from nucleus (in Bohr radii)
    Z : int
        Nuclear charge (Z=1 for hydrogen)
    
    Returns:
    --------
    R_nl : ndarray
        Radial wavefunction values
    """
    # Normalization constant
    rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n)**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    
    # Associated Laguerre polynomial
    laguerre = genlaguerre(n - l - 1, 2 * l + 1)
    
    # Radial wavefunction
    R_nl = norm * np.exp(-rho / 2) * rho**l * laguerre(rho)
    
    return R_nl

def angular_wavefunction(l, m, theta, phi):
    """
    Calculate the angular wavefunction Y_lm(θ,φ) (spherical harmonic).
    
    Parameters:
    -----------
    l : int
        Angular momentum quantum number
    m : int
        Magnetic quantum number (-l ≤ m ≤ l)
    theta : ndarray
        Polar angle (0 to π)
    phi : ndarray
        Azimuthal angle (0 to 2π)
    
    Returns:
    --------
    Y_lm : ndarray (complex)
        Angular wavefunction (spherical harmonic)
    """
    # scipy.special.sph_harm uses (m, l, phi, theta) convention
    return sph_harm(m, l, phi, theta)

def total_wavefunction(n, l, m, r, theta, phi, Z=1):
    """
    Calculate the total wavefunction Ψ_nlm(r,θ,φ) = R_nl(r) * Y_lm(θ,φ).
    
    Returns:
    --------
    psi : ndarray (complex)
        Total wavefunction
    """
    R_nl = radial_wavefunction(n, l, r, Z)
    Y_lm = angular_wavefunction(l, m, theta, phi)
    return R_nl * Y_lm

def probability_density(n, l, m, r, theta, phi, Z=1):
    """
    Calculate probability density |Ψ|² for orbital visualization.
    """
    psi = total_wavefunction(n, l, m, r, theta, phi, Z)
    return np.abs(psi)**2

def orbital_energy(n, Z=1):
    """
    Calculate orbital energy for hydrogen-like atoms.
    
    E_n = -Z²·13.6 eV / n²
    
    Parameters:
    -----------
    n : int
        Principal quantum number
    Z : int
        Nuclear charge
    
    Returns:
    --------
    energy : float
        Orbital energy in eV
    """
    return -Z**2 * RYDBERG_ENERGY / n**2

def orbital_name(n, l, m):
    """
    Convert quantum numbers to orbital name (e.g., 1s, 2p_x, 3d_z²).
    """
    subshell_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    
    if l not in subshell_names:
        return f"{n}(l={l})"
    
    base_name = f"{n}{subshell_names[l]}"
    
    # Add subscript for p, d, f orbitals
    if l == 1:  # p orbitals
        p_names = {-1: 'y', 0: 'z', 1: 'x'}
        return f"{base_name}_{p_names.get(m, m)}"
    elif l == 2:  # d orbitals
        d_names = {-2: 'xy', -1: 'yz', 0: 'z²', 1: 'xz', 2: 'x²-y²'}
        return f"{base_name}_{d_names.get(m, m)}"
    elif l == 3:  # f orbitals
        return f"{base_name}_m={m}"
    
    return base_name

def create_spherical_grid(n_theta=100, n_phi=100, r_val=1.0):
    """
    Create a spherical coordinate grid for 3D plotting.
    """
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta, phi = np.meshgrid(theta, phi)
    
    return theta, phi

def spherical_to_cartesian(r, theta, phi):
    """
    Convert spherical coordinates to Cartesian.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

# ============================================================================
# MAIN DEMONSTRATION
# ============================================================================

print("=" * 80)
print("ATOMIC ORBITALS: SHAPES AND ENERGIES")
print("=" * 80)
print("\nDemonstrating the quantum mechanical description of atomic orbitals")
print("for hydrogen-like atoms.\n")

# Create comprehensive figure with multiple subplots
fig = plt.figure(figsize=(20, 14))
gs = GridSpec(4, 5, figure=fig, hspace=0.35, wspace=0.35)

# ============================================================================
# PART 1: ENERGY LEVEL DIAGRAM
# ============================================================================

ax1 = fig.add_subplot(gs[0, 0])

n_levels = [1, 2, 3, 4]
energies = [orbital_energy(n) for n in n_levels]

# Plot energy levels
y_offset = 0
for i, n in enumerate(n_levels):
    E = energies[i]
    
    # Draw energy level line
    max_orbitals = n  # Number of subshells for this n
    line_width = max_orbitals * 0.3
    
    ax1.plot([0.5 - line_width/2, 0.5 + line_width/2], [E, E], 
            'k-', linewidth=2)
    
    # Label the level
    ax1.text(0.5 + line_width/2 + 0.1, E, f'n={n}\n{E:.2f} eV', 
            va='center', fontsize=9)
    
    # Show subshells
    subshell_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    x_pos = 0.5 - line_width/2
    for l in range(n):
        x_pos += 0.15
        subshell = subshell_names.get(l, f'l={l}')
        ax1.text(x_pos, E - 0.3, f'{subshell}', fontsize=8, ha='center')

ax1.set_xlim(0, 2)
ax1.set_ylim(min(energies) - 1, 0.5)
ax1.set_ylabel('Energy (eV)', fontsize=10, fontweight='bold')
ax1.set_title('Hydrogen Atom Energy Levels', fontsize=11, fontweight='bold')
ax1.axhline(y=0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='Ionization (E=0)')
ax1.grid(True, alpha=0.3, axis='y')
ax1.set_xticks([])
ax1.legend(fontsize=8)

# ============================================================================
# PART 2: RADIAL PROBABILITY DISTRIBUTIONS
# ============================================================================

ax2 = fig.add_subplot(gs[0, 1:3])

r = np.linspace(0.01, 20, 500)

orbitals_to_plot = [
    (1, 0, '1s', '#FF6B6B'),
    (2, 0, '2s', '#4ECDC4'),
    (2, 1, '2p', '#45B7D1'),
    (3, 0, '3s', '#FFA07A'),
    (3, 1, '3p', '#98D8C8'),
    (3, 2, '3d', '#95E1D3'),
]

for n, l, name, color in orbitals_to_plot:
    R = radial_wavefunction(n, l, r)
    # Radial probability density: 4πr²|R_nl|²
    prob_density = 4 * np.pi * r**2 * R**2
    ax2.plot(r, prob_density, label=name, linewidth=2, color=color)

ax2.set_xlabel('Distance from nucleus (Bohr radii)', fontsize=10, fontweight='bold')
ax2.set_ylabel('Radial probability density', fontsize=10, fontweight='bold')
ax2.set_title('Radial Probability Distributions', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9, ncol=2)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 20)

# ============================================================================
# PART 3: s ORBITAL (n=1, l=0)
# ============================================================================

ax3 = fig.add_subplot(gs[0, 3], projection='3d')

# Create spherical grid
theta, phi = create_spherical_grid(50, 50)

# For s orbital, angular part is constant (spherically symmetric)
# Use radial probability to determine surface radius
r_surface = 1.5  # Bohr radii
prob_s = probability_density(1, 0, 0, r_surface, theta, phi)

# Convert to Cartesian
x, y, z = spherical_to_cartesian(r_surface, theta, phi)

# Color by probability
colors_s = prob_s / prob_s.max()

surf = ax3.plot_surface(x, y, z, facecolors=cm.Reds(colors_s), 
                        alpha=0.8, edgecolor='none')

ax3.set_xlabel('x (a₀)', fontsize=8)
ax3.set_ylabel('y (a₀)', fontsize=8)
ax3.set_zlabel('z (a₀)', fontsize=8)
ax3.set_title('1s Orbital\n(Spherical)', fontsize=10, fontweight='bold')
ax3.set_box_aspect([1, 1, 1])
ax3.view_init(elev=20, azim=45)

# ============================================================================
# PART 4: 2s ORBITAL (showing radial nodes)
# ============================================================================

ax4 = fig.add_subplot(gs[0, 4], projection='3d')

r_surface_2s = 5.0
prob_2s = probability_density(2, 0, 0, r_surface_2s, theta, phi)

x_2s, y_2s, z_2s = spherical_to_cartesian(r_surface_2s, theta, phi)
colors_2s = prob_2s / prob_2s.max()

ax4.plot_surface(x_2s, y_2s, z_2s, facecolors=cm.Blues(colors_2s), 
                alpha=0.7, edgecolor='none')

# Add inner node sphere
r_node = 2.0
x_node, y_node, z_node = spherical_to_cartesian(r_node, theta, phi)
ax4.plot_surface(x_node, y_node, z_node, color='lightblue', 
                alpha=0.3, edgecolor='k', linewidth=0.3)

ax4.set_xlabel('x (a₀)', fontsize=8)
ax4.set_ylabel('y (a₀)', fontsize=8)
ax4.set_zlabel('z (a₀)', fontsize=8)
ax4.set_title('2s Orbital\n(1 Radial Node)', fontsize=10, fontweight='bold')
ax4.set_box_aspect([1, 1, 1])
ax4.view_init(elev=20, azim=45)

# ============================================================================
# PART 5: p ORBITALS (l=1, m=-1,0,1) - 3D surfaces
# ============================================================================

p_orbitals = [
    (2, 1, 1, gs[1, 0], 'p_x', '#4ECDC4'),
    (2, 1, 0, gs[1, 1], 'p_z', '#45B7D1'),
    (2, 1, -1, gs[1, 2], 'p_y', '#95E1D3'),
]

for n, l, m, grid_pos, name, color in p_orbitals:
    ax = fig.add_subplot(grid_pos, projection='3d')
    
    # Create probability density
    r_p = 4.0
    prob_p = np.real(probability_density(n, l, m, r_p, theta, phi))
    
    # Scale radius by probability for visualization
    r_plot = r_p * (1 + 0.3 * prob_p / prob_p.max())
    
    x_p, y_p, z_p = spherical_to_cartesian(r_plot, theta, phi)
    
    # Color by sign of wavefunction
    psi_p = np.real(total_wavefunction(n, l, m, r_p, theta, phi))
    colors = np.where(psi_p > 0, 0.8, 0.2)
    
    ax.plot_surface(x_p, y_p, z_p, facecolors=cm.cool(colors), 
                   alpha=0.85, edgecolor='none')
    
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.set_zlabel('z', fontsize=8)
    ax.set_title(f'2{name} Orbital\n(Dumbbell)', fontsize=10, fontweight='bold')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=20, azim=45)

# ============================================================================
# PART 6: d ORBITALS (l=2) - 5 different shapes
# ============================================================================

d_orbitals = [
    (3, 2, 0, gs[1, 3], 'd_z²', '#FFA07A'),
    (3, 2, 2, gs[1, 4], 'd_x²-y²', '#FF7F50'),
    (3, 2, 1, gs[2, 0], 'd_xz', '#FF6347'),
    (3, 2, -1, gs[2, 1], 'd_yz', '#FFB6C1'),
    (3, 2, -2, gs[2, 2], 'd_xy', '#FFC0CB'),
]

for n, l, m, grid_pos, name, color in d_orbitals:
    ax = fig.add_subplot(grid_pos, projection='3d')
    
    r_d = 6.0
    prob_d = np.real(probability_density(n, l, m, r_d, theta, phi))
    
    # Scale radius
    r_plot = r_d * (1 + 0.4 * prob_d / (prob_d.max() + 1e-10))
    
    x_d, y_d, z_d = spherical_to_cartesian(r_plot, theta, phi)
    
    # Color by wavefunction sign
    psi_d = np.real(total_wavefunction(n, l, m, r_d, theta, phi))
    colors_d = np.where(psi_d > 0, 0.8, 0.2)
    
    ax.plot_surface(x_d, y_d, z_d, facecolors=cm.autumn(colors_d), 
                   alpha=0.8, edgecolor='none')
    
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.set_zlabel('z', fontsize=8)
    ax.set_title(f'3{name} Orbital\n(Cloverleaf)', fontsize=10, fontweight='bold')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=20, azim=45)

# ============================================================================
# PART 7: f ORBITALS (l=3) - Examples
# ============================================================================

f_orbitals = [
    (4, 3, 0, gs[2, 3], 'f (m=0)', '#DDA0DD'),
    (4, 3, 3, gs[2, 4], 'f (m=3)', '#DA70D6'),
    (4, 3, -3, gs[3, 0], 'f (m=-3)', '#BA55D3'),
]

for n, l, m, grid_pos, name, color in f_orbitals:
    ax = fig.add_subplot(grid_pos, projection='3d')
    
    r_f = 8.0
    prob_f = np.real(probability_density(n, l, m, r_f, theta, phi))
    
    r_plot = r_f * (1 + 0.5 * prob_f / (prob_f.max() + 1e-10))
    
    x_f, y_f, z_f = spherical_to_cartesian(r_plot, theta, phi)
    
    psi_f = np.real(total_wavefunction(n, l, m, r_f, theta, phi))
    colors_f = np.where(psi_f > 0, 0.8, 0.2)
    
    ax.plot_surface(x_f, y_f, z_f, facecolors=cm.Purples(colors_f), 
                   alpha=0.75, edgecolor='none')
    
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.set_zlabel('z', fontsize=8)
    ax.set_title(f'4{name} Orbital\n(Complex)', fontsize=10, fontweight='bold')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=20, azim=45)

# ============================================================================
# PART 8: CONTOUR PLOTS (Cross-sections)
# ============================================================================

# 2p_z contour in x-z plane
ax8 = fig.add_subplot(gs[3, 1])

x_range = np.linspace(-8, 8, 200)
z_range = np.linspace(-8, 8, 200)
X_xz, Z_xz = np.meshgrid(x_range, z_range)

# Convert to spherical coordinates (y=0 plane)
R_xz = np.sqrt(X_xz**2 + Z_xz**2)
Theta_xz = np.arctan2(np.sqrt(X_xz**2), Z_xz)
Phi_xz = np.arctan2(0, X_xz)

# Calculate wavefunction
psi_2pz = np.real(total_wavefunction(2, 1, 0, R_xz, Theta_xz, Phi_xz))

levels = np.linspace(-0.15, 0.15, 20)
contour = ax8.contourf(X_xz, Z_xz, psi_2pz, levels=levels, cmap='RdBu_r')
ax8.contour(X_xz, Z_xz, psi_2pz, levels=[0], colors='black', linewidths=2)

ax8.set_xlabel('x (a₀)', fontsize=9, fontweight='bold')
ax8.set_ylabel('z (a₀)', fontsize=9, fontweight='bold')
ax8.set_title('2p_z Contour (x-z plane)\nBlue: ψ>0, Red: ψ<0', fontsize=10, fontweight='bold')
ax8.set_aspect('equal')
plt.colorbar(contour, ax=ax8, label='ψ')

# 3d_z² contour in x-z plane
ax9 = fig.add_subplot(gs[3, 2])

psi_3dz2 = np.real(total_wavefunction(3, 2, 0, R_xz, Theta_xz, Phi_xz))

levels_d = np.linspace(-0.05, 0.05, 20)
contour_d = ax9.contourf(X_xz, Z_xz, psi_3dz2, levels=levels_d, cmap='RdBu_r')
ax9.contour(X_xz, Z_xz, psi_3dz2, levels=[0], colors='black', linewidths=2)

ax9.set_xlabel('x (a₀)', fontsize=9, fontweight='bold')
ax9.set_ylabel('z (a₀)', fontsize=9, fontweight='bold')
ax9.set_title('3d_z² Contour (x-z plane)\nShowing Nodal Surfaces', fontsize=10, fontweight='bold')
ax9.set_aspect('equal')
plt.colorbar(contour_d, ax=ax9, label='ψ')

# ============================================================================
# PART 9: ANGULAR DISTRIBUTION PATTERNS
# ============================================================================

ax10 = fig.add_subplot(gs[3, 3], projection='3d')

# Plot angular distribution for d_z² orbital
theta_ang = np.linspace(0, np.pi, 100)
phi_ang = np.linspace(0, 2 * np.pi, 100)
theta_ang, phi_ang = np.meshgrid(theta_ang, phi_ang)

# Get angular part only
Y_lm_d = angular_wavefunction(2, 0, theta_ang, phi_ang)
r_ang = np.abs(Y_lm_d)**2 * 10  # Scale for visibility

x_ang, y_ang, z_ang = spherical_to_cartesian(r_ang, theta_ang, phi_ang)

# Color by sign
colors_ang = np.where(np.real(Y_lm_d) > 0, 0.8, 0.2)

ax10.plot_surface(x_ang, y_ang, z_ang, facecolors=cm.RdYlBu(colors_ang), 
                 alpha=0.8, edgecolor='none')

ax10.set_xlabel('x', fontsize=8)
ax10.set_ylabel('y', fontsize=8)
ax10.set_zlabel('z', fontsize=8)
ax10.set_title('Angular Distribution\nd_z² (|Y₂₀|²)', fontsize=10, fontweight='bold')
ax10.set_box_aspect([1, 1, 1])
ax10.view_init(elev=20, azim=45)

# ============================================================================
# PART 10: KEY PRINCIPLES TEXT
# ============================================================================

ax11 = fig.add_subplot(gs[3, 4])

principles_text = """
ORBITAL CHARACTERISTICS

s orbitals (l=0):
• Spherically symmetric
• No angular nodes
• n-1 radial nodes

p orbitals (l=1):
• Dumbbell shaped
• 1 angular node (plane)
• n-2 radial nodes
• 3 orientations (px, py, pz)

d orbitals (l=2):
• Cloverleaf patterns
• 2 angular nodes
• n-3 radial nodes
• 5 orientations

f orbitals (l=3):
• Complex shapes
• 3 angular nodes
• n-4 radial nodes
• 7 orientations

NODE RULES:
Total nodes = n - 1
Angular nodes = l
Radial nodes = n - l - 1

ENERGY (H-atom):
E_n = -13.6 eV / n²
(Energy depends only on n,
not on l or m)
"""

ax11.text(0.05, 0.95, principles_text, 
         ha='left', va='top', fontsize=8.5,
         transform=ax11.transAxes, family='monospace')
ax11.axis('off')

# ============================================================================
# ADDITIONAL ANALYSIS PLOT
# ============================================================================

ax12 = fig.add_subplot(gs[3, 0])

# Show number of orbitals per shell
shells = [1, 2, 3, 4]
s_orbitals = [1, 1, 1, 1]
p_orbitals = [0, 3, 3, 3]
d_orbitals = [0, 0, 5, 5]
f_orbitals = [0, 0, 0, 7]

x = np.arange(len(shells))
width = 0.6

bottom1 = np.array(s_orbitals)
bottom2 = bottom1 + np.array(p_orbitals)
bottom3 = bottom2 + np.array(d_orbitals)

ax12.bar(x, s_orbitals, width, label='s (l=0)', color='#FF6B6B')
ax12.bar(x, p_orbitals, width, bottom=bottom1, label='p (l=1)', color='#4ECDC4')
ax12.bar(x, d_orbitals, width, bottom=bottom2, label='d (l=2)', color='#FFA07A')
ax12.bar(x, f_orbitals, width, bottom=bottom3, label='f (l=3)', color='#DDA0DD')

ax12.set_xlabel('Shell (n)', fontsize=9, fontweight='bold')
ax12.set_ylabel('Number of Orbitals', fontsize=9, fontweight='bold')
ax12.set_title('Orbitals per Shell\nTotal = n²', fontsize=10, fontweight='bold')
ax12.set_xticks(x)
ax12.set_xticklabels([f'n={n}' for n in shells])
ax12.legend(fontsize=8)
ax12.grid(True, alpha=0.3, axis='y')

# Add total count labels
totals = [s + p + d + f for s, p, d, f in zip(s_orbitals, p_orbitals, d_orbitals, f_orbitals)]
for i, total in enumerate(totals):
    ax12.text(i, total + 0.5, f'{total}', ha='center', fontweight='bold', fontsize=9)

plt.suptitle('Atomic Orbitals: Shapes and Energies (s, p, d, f)', 
            fontsize=18, fontweight='bold', y=0.998)

plt.show(block=True)

# ============================================================================
# TERMINAL OUTPUT SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("SUMMARY: ATOMIC ORBITAL CHARACTERISTICS")
print("=" * 80)

orbital_data = [
    ('1s', 1, 0, orbital_energy(1), 0, 0),
    ('2s', 2, 0, orbital_energy(2), 1, 0),
    ('2p', 2, 1, orbital_energy(2), 0, 1),
    ('3s', 3, 0, orbital_energy(3), 2, 0),
    ('3p', 3, 1, orbital_energy(3), 1, 1),
    ('3d', 3, 2, orbital_energy(3), 0, 2),
    ('4s', 4, 0, orbital_energy(4), 3, 0),
    ('4p', 4, 1, orbital_energy(4), 2, 1),
    ('4d', 4, 2, orbital_energy(4), 1, 2),
    ('4f', 4, 3, orbital_energy(4), 0, 3),
]

print(f"\n{'Orbital':<8} {'n':<3} {'l':<3} {'Energy (eV)':<15} {'Radial Nodes':<15} {'Angular Nodes':<15}")
print("-" * 80)

for name, n, l, energy, rad_nodes, ang_nodes in orbital_data:
    print(f"{name:<8} {n:<3} {l:<3} {energy:>12.3f}   {rad_nodes:^15} {ang_nodes:^15}")

print("\n" + "=" * 80)
print("KEY INSIGHTS:")
print("=" * 80)
print("• Orbitals are solutions to the Schrödinger equation: Ψ(r,θ,φ) = R_nl(r) Y_lm(θ,φ)")
print("• Energy depends only on n for hydrogen: E_n = -13.6 eV / n²")
print("• Orbital shape determined by l (angular momentum quantum number):")
print("    - l=0 (s): Spherical")
print("    - l=1 (p): Dumbbell (3 orientations)")
print("    - l=2 (d): Cloverleaf (5 orientations)")
print("    - l=3 (f): Complex (7 orientations)")
print("• Number of orbitals in shell n: n² total orbitals")
print("• Nodes (zero-probability regions):")
print("    - Total nodes = n - 1")
print("    - Angular nodes (nodal planes/cones) = l")
print("    - Radial nodes (nodal spheres) = n - l - 1")
print("• Orbitals with same n but different l are degenerate in H-atom")
print("• In multi-electron atoms, energy depends on both n and l")
print("=" * 80)
print("\nVisualization complete! Interactive plot displayed.")
print("=" * 80)

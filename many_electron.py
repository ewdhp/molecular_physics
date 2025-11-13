import numpy as np
import matplotlib.pyplot as plt

# Slater's rules for effective nuclear charge (Z_eff)
def slater_zeff(Z, n, electrons):
    """
    Calculate effective nuclear charge for a given electron using Slater's rules.
    Z: atomic number
    n: principal quantum number of electron of interest
    electrons: list of (n, count) for all electrons in atom
    """
    S = 0
    for n_e, count in electrons:
        if n_e < n:
            S += 0.85 * count  # electrons in shells below
        elif n_e == n:
            S += 0.35 * (count - 1)  # other electrons in same shell
    return Z - S

# Example: Calculate Z_eff for 3p electron in phosphorus (Z=15)
Z = 15
electrons = [(1, 2), (2, 8), (3, 5)]  # (n, count)
n_interest = 3
zeff = slater_zeff(Z, n_interest, electrons)
print(f"Effective nuclear charge (Z_eff) for 3p electron in phosphorus: {zeff:.2f}")

# Plot Z_eff for 2p electrons across period 2 elements (Li to Ne)
Z_list = np.arange(3, 11)
zeff_list = []
for Z in Z_list:
    # Electron configuration for period 2: 1s^2 2s^2 2p^(Z-4)
    electrons = [(1, 2), (2, Z-2)]
    zeff = slater_zeff(Z, 2, electrons)
    zeff_list.append(zeff)

elements = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
plt.figure(figsize=(8,5))
plt.plot(elements, zeff_list, marker='o')
plt.title('Effective Nuclear Charge (Z_eff) for 2p Electrons Across Period 2')
plt.xlabel('Element')
plt.ylabel('Z_eff')
plt.grid(True)
plt.show()
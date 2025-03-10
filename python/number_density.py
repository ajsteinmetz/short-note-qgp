import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kn, zeta

# Constants
g_mu = 2
m_mu_c2 = 105.66  # Muon rest mass energy in MeV
hbar_c = 197.327  # MeV*fm

# Conversion factor (1 fm = 1e-13 cm)
fm_to_cm = 1e-13

# Functions for number densities
def n_mu_1(T):
    factor = g_mu / (2 * np.pi**2)
    x = m_mu_c2 / T
    prefactor = (T / (hbar_c))**3
    return factor * prefactor * x**2 * kn(2, x) / fm_to_cm**3

def n_mu_2(T):
    factor = g_mu
    prefactor = ((m_mu_c2 * T) / (2 * np.pi * (hbar_c)**2))**(3/2)
    exponent = np.exp(-m_mu_c2 / T)
    return factor * prefactor * exponent / fm_to_cm**3

def n_mu_3(T):
    factor = (1.10921) * g_mu * (3 / 4) * zeta(3) / (np.pi**2)
    prefactor = (T / (hbar_c))**3
    return factor * prefactor / fm_to_cm**3

# Temperature range in MeV
T = np.linspace(145000, 150000, 1000)

# Calculate number densities
n1 = n_mu_1(T)
n2 = n_mu_2(T)
n3 = n_mu_3(T)

# Plotting number densities
plt.figure(figsize=(10,6))
plt.plot(T, n1, label=r'$n^{1}_{\mu}$', linewidth=2)
plt.plot(T, n2, '--', label=r'$n^{2}_{\mu}$', linewidth=2)
plt.plot(T, n3, '-.', label=r'$n^{3}_{\mu}$', linewidth=2)
plt.xlabel(r'$k_{B}T\,(\mathrm{MeV})$', fontsize=14)
plt.ylabel(r'Muon number density $(\mathrm{cm}^{-3})$', fontsize=14)
plt.title('Muon Number Density vs Temperature', fontsize=16)
plt.legend()
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.show()

# Plotting relative differences between number densities
plt.figure(figsize=(10,6))
# plt.plot(T, np.abs(n2 - n1)/n1, label=r'$|(n^{1}_{\mu}-n^{2}_{\mu})/n^{1}_{\mu}|$', color='purple', linewidth=2)
plt.plot(T, np.abs(n3 - n1)/n1, label=r'$|(n^{1}_{\mu}-n^{3}_{\mu})/n^{1}_{\mu}|$', color='green', linewidth=2, linestyle='--')
plt.xlabel(r'$k_{B}T\,(\mathrm{MeV})$', fontsize=14)
plt.ylabel(r'Relative Difference', fontsize=14)
plt.title('Relative Differences in Muon Number Densities vs Temperature', fontsize=16)
plt.legend()
plt.grid()
plt.yscale('log')
# plt.xscale('log')
plt.tight_layout()
plt.show()

# Print relative difference at the final temperature
final_relative_difference_n1_n3 = np.abs(n1[-1] - n3[-1]) / n1[-1]
print(f'Relative difference at final temperature: {final_relative_difference_n1_n3:.4e}')

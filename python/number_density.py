import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.special import kv  # Modified Bessel function of the second kind

class particle:
    def __init__(self, name, mass, charge, degeneracy):
        """
        Initialize a particle instance.
        
        Parameters:
            name (str): Name of the particle.
            mass (float): Mass in MeV.
            charge (float): Electric charge (in units of e).
            degeneracy (int): Degrees of freedom.
        """
        self.name = name
        self.mass = mass        # in MeV
        self.charge = charge    # in units of e
        self.degeneracy = degeneracy

    @property
    def magneton(self):
        """
        Returns the magnetic moment in units of the Bohr magneton.
        
        Assumes a gyromagnetic factor g=2 so that:
            μ/μ_B = (g/2) * (q * m_e / m) = q * m_e / m.
        The electron mass m_e is taken as 0.511 MeV.
        """
        m_electron = 0.511  # in MeV
        return abs(self.charge) * m_electron / self.mass

    def __repr__(self):
        return (f"particle(name='{self.name}', mass={self.mass} MeV, "
                f"charge={self.charge}e, degeneracy={self.degeneracy}, "
                f"magneton={self.magneton:.6f} μ_B)")

def boltzmann_number_density(m, T, g, num_terms=10):
    """
    Compute the number density for a heavy particle using the Boltzmann expansion 
    of the Fermi-Dirac distribution.
    
    The expression is:
        n(T) = (g/(2π²)) * T³ * ∑ₖ₌₁^(num_terms) [ (-1)^(k+1) / k⁴ * (k*m/T)² * K₂(k*m/T) ]
    Units are chosen so that k_B = ħ = c = 1.
    """
    prefactor = g / (2 * np.pi**2) * T**3
    series_sum = 0.0
    for k in range(1, num_terms + 1):
        arg = k * m / T
        term = ((-1)**(k+1) / k**4) * (k * m / T)**2 * kv(2, arg)
        series_sum += term
    return prefactor * series_sum

def massless_number_density(T, g):
    """
    Compute the number density for a massless fermion gas.
    
    n(T) = (g/(2π²)) * (3ζ(3)/2) * T³,
    with ζ(3) ≃ 1.202.
    """
    zeta3 = 1.202
    return g / (2 * np.pi**2) * (3 * zeta3 / 2) * T**3

def heavy_quark_number_density(m, T, g, num_terms=10):
    """
    Compute the number density for quarks with a threshold interpolation.
    
    For T ≤ 150 MeV: n(T) = 0.
    For 150 MeV < T < 170 MeV: 
         n(T) = ((T - 150) / 20) * n_boltz(T),
         where n_boltz(T) is the Boltzmann number density.
    For T ≥ 170 MeV: n(T) = n_boltz(T).
    """
    T = np.array(T)
    n = np.empty_like(T)
    mask_low = T <= 150
    mask_trans = (T > 150) & (T < 170)
    mask_high = T >= 170

    n[mask_low] = 0.0
    n[mask_trans] = ((T[mask_trans] - 150) / 20) * boltzmann_number_density(m, T[mask_trans], g, num_terms)
    n[mask_high] = boltzmann_number_density(m, T[mask_high], g, num_terms)
    return n

# Constants and conversion factors
T0 = 2.35e-4             # Today's CMB temperature in eV
B0_low = 1e-20           # Lower end of today's intergalactic magnetic field (Tesla)
B0_high = 1e-12          # Upper end (Tesla)
B_critical = 4.41e9      # Schwinger critical field in Tesla
MeV2_to_Tesla = 5.122e9  # Conversion factor from (MeV)² to Tesla
mu_B = 0.296             # Bohr magneton in MeV⁻¹

# Define particles using the provided information.
# All particles use the heavy (Boltzmann) number density.
# Leptons:
electron = particle("electron", 0.511, -1, 2)
muon     = particle("muon", 105.7, -1, 2)
tau      = particle("tau", 1776.9, -1, 2)

# Quarks: all quarks use the heavy_quark_number_density with interpolation.
up       = particle("up", 2.3, 2/3, 6)
down     = particle("down", 4.8, -1/3, 6)
strange  = particle("strange", 96.0, -1/3, 6)
charm    = particle("charm", 1270.0, 2/3, 6)
bottom   = particle("bottom", 4180.0, -1/3, 6)

# Temperature range in MeV; 1000 points between 5 and 500 MeV.
T_lo = 5
T_hi = 500
T_MeV = np.linspace(T_lo, T_hi, 1000)
T_eV = T_MeV * 1e6  # Convert temperature from MeV to eV

# Compute the primordial magnetic field using the scaling law:
# B_PMF(T) = B0 * (T²)/(T0²)
B_low = B0_low * (T_eV**2) / (T0**2)
B_high = B0_high * (T_eV**2) / (T0**2)

# --- Compute magnetizations ---
# Magnetization is defined as: 
#   B_i(T) = μ_B * (μ_i/μ_B) * n_i(T) = μ_i * n_i(T)
# (The μ_B factor converts units; here μ_B = 0.296 MeV⁻¹.)

# Leptons: use heavy number density for all.
n_electron = boltzmann_number_density(electron.mass, T_MeV, electron.degeneracy)
n_muon = boltzmann_number_density(muon.mass, T_MeV, muon.degeneracy)
n_tau  = boltzmann_number_density(tau.mass, T_MeV, tau.degeneracy)

B_electron = 2 * MeV2_to_Tesla * mu_B * electron.magneton * n_electron
B_muon = 2 * MeV2_to_Tesla * mu_B * muon.magneton * n_muon
B_tau  = 2 * MeV2_to_Tesla * mu_B * tau.magneton * n_tau

# Quarks: use heavy_quark_number_density with interpolation.
n_up = heavy_quark_number_density(up.mass, T_MeV, up.degeneracy)
n_down = heavy_quark_number_density(down.mass, T_MeV, down.degeneracy)
n_strange = heavy_quark_number_density(strange.mass, T_MeV, strange.degeneracy)
n_charm = heavy_quark_number_density(charm.mass, T_MeV, charm.degeneracy)
n_bottom = heavy_quark_number_density(bottom.mass, T_MeV, bottom.degeneracy)

B_up = 2 * MeV2_to_Tesla * mu_B * up.magneton * n_up
B_down = 2 * MeV2_to_Tesla * mu_B * down.magneton * n_down
B_strange = 2 * MeV2_to_Tesla * mu_B * strange.magneton * n_strange
B_charm = 2 * MeV2_to_Tesla * mu_B * charm.magneton * n_charm
B_bottom = 2 * MeV2_to_Tesla * mu_B * bottom.magneton * n_bottom

# --- Plotting ---
fig, ax = plt.subplots(figsize=(10, 6))

# Shade the region between the lower and upper bounds for the primordial magnetic field.
ax.fill_between(T_MeV, B_low, B_high, color='lightgrey', alpha=1.0, label='Allowed PMF range')
# Add a vertical line at T = 150 MeV to indicate the start of the transition region,
ax.axvline(x=150, color='green', linestyle='-', linewidth=1.0, alpha=1.0, label=r'$T=150\,\mathrm{MeV}$')

# Color grading for leptons (using different shades of blue).
colors_leptons = {
    'electron': plt.cm.Blues(1.00),
    'muon': plt.cm.Blues(0.75),
    'tau': plt.cm.Blues(0.50)
}

ax.plot(T_MeV, B_electron, label='Electron Magnetization', color=colors_leptons['electron'], linestyle='solid')
ax.plot(T_MeV, B_muon, label='Muon Magnetization', color=colors_leptons['muon'], linestyle='solid')
ax.plot(T_MeV, B_tau, label='Tau Magnetization', color=colors_leptons['tau'], linestyle='solid')

# Color grading for quarks (using different shades of red).
colors_quarks = {
    'light': plt.cm.Reds(0.4),
    'strange': plt.cm.Reds(0.6),
    'charm': plt.cm.Reds(0.8),
    'bottom': plt.cm.Reds(1.0)
}

# For quarks, split the curve into transition (dashed) and non-transition (solid) segments.
mask_trans = (T_MeV > 150) & (T_MeV < 170)
mask_non_trans = (T_MeV >= 170)

def plot_quark(T, B, label, color):
    # Plot non-transition segments with a solid line (include label only once).
    ax.plot(T[mask_non_trans], B[mask_non_trans], label=label, color=color, linestyle='solid')
    # Plot the transition region with a dashed line (no label to avoid duplicate legend entry).
    ax.plot(T[mask_trans], B[mask_trans], color=color, linestyle='dashed')

plot_quark(T_MeV, B_up + B_down, 'Light-Quark Magnetization', colors_quarks['light'])
plot_quark(T_MeV, B_strange, 'Strange-Quark Magnetization', colors_quarks['strange'])
plot_quark(T_MeV, B_charm, 'Charm-Quark Magnetization', colors_quarks['charm'])
plot_quark(T_MeV, B_bottom, 'Bottom-Quark Magnetization', colors_quarks['bottom'])

ax.set_xlabel("Temperature (MeV)", fontsize=14)
ax.set_ylabel("Magnetic Field (T)", fontsize=14)
ax.set_yscale('log')
ax.set_xscale('log')
ax.invert_xaxis()  # Universe gets cooler from left to right
ax.set_xlim(T_hi, T_lo)
ax.set_ylim(1e0, 1e18)  # y-axis from 10^0 to 10^18 Tesla

# Major gridlines (already dashed) for both axes.
ax.grid(which='major', linestyle='--', linewidth=0.75, alpha=1.0)
# Add dashed gridlines for the minor x-axis.
ax.grid(which='minor', axis='x', linestyle='--', linewidth=0.75, alpha=1.0)

ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=100))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=100))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())

#ax.legend()

plt.show()

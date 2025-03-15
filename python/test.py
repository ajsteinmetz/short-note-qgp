#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def fermi(x):
    """Fermi-Dirac distribution function."""
    return 1.0 / (np.exp(x) + 1.0)

def integrand(E, s, Omega, T, m, B, e, g):
    """
    Computes the integrand for the dimensionless energy variable E,
    for a given spin state s, spin potential Omega, temperature T,
    mass m, magnetic field B, electron charge e, and g-factor g.
    """
    # Effective mass squared: m^2(s) = m^2 - e*B*g*s.
    m_eff_sq = m**2 - e * B * g * s
    if m_eff_sq < 0:
        # Avoid negative square roots if B is too high for this spin state.
        return 0.0
    # Dimensionless effective mass: m_T = sqrt(m_eff_sq)/T.
    m_T = np.sqrt(m_eff_sq) / T
    # (E^2 - m_T^2)^(3/2)
    factor = (E**2 - m_T**2)**(3/2)
    
    # With eta = 0, the effective dimensionless chemical potential is:
    # Sigma_T = s*Omega/T.
    Sigma_T = s * Omega / T
    F_val = fermi(E - Sigma_T)
    
    # Return the integrand as given in eq. (mag_form).
    return factor * (F_val * (1 - F_val) / E + F_val / E**2)

def compute_integral(s, Omega, T, m, B, e, g):
    """
    Computes the integral I(s) for a given spin state s,
    from E = m_T(s) to infinity.
    """
    m_eff_sq = m**2 - e * B * g * s
    if m_eff_sq < 0:
        return 0.0
    m_T = np.sqrt(m_eff_sq) / T
    result, error = quad(integrand, m_T, np.inf, args=(s, Omega, T, m, B, e, g), limit=100)
    return result

def compute_magnetization(T, m, g, e, B, Omega, n_C=1.0, V=1.0):
    """
    Computes the magnetization M at temperature T.
    
    Parameters:
      T    : Temperature in MeV.
      m    : Electron mass in MeV.
      g    : g-factor.
      e    : Absolute electron charge in natural units.
      B    : Magnetic field strength.
      Omega: Spin polarization potential in MeV.
      n_C  : Color factor (electrons are colorless, so n_C=1).
      V    : Volume (set to 1 in arbitrary units).
      
    Returns:
      M    : Magnetization in natural units.
    """
    # Compute the integral for each spin state.
    spins = [0.5, -0.5]
    I_vals = {s: compute_integral(s, Omega, T, m, B, e, g) for s in spins}
    # The overall spin sum gives I(0.5) - I(-0.5).
    I_total = I_vals[0.5] - I_vals[-0.5]
    
    # Overall prefactor from eq. (mag_form):
    #   -(g|Q|)/(2T) * (2 n_C V T^3)/(3(2π)^2)
    prefactor = - (g * e) / (2 * T) * (2 * n_C * V * T**3) / (3 * (2*np.pi)**2)
    M = prefactor * I_total
    return M

def main():
    # Electron and simulation parameters (in natural units)
    m = 0.511      # Electron mass in MeV
    g = 2.0        # g-factor for the electron
    e_val = 0.303  # Absolute electron charge in natural units
    B = 0.01       # Magnetic field strength in natural units (e.g., MeV^2)
    Omega = 0.0   # Spin polarization potential in MeV
    MeV2_to_Tesla = 5.122e9  # Conversion factor from (MeV)² to Tesla

    # Define the temperature range: from 5 MeV to 500 MeV.
    T_vals = np.linspace(5, 500, 100)
    M_vals = []

    # Compute magnetization for each temperature.
    for T in T_vals:
        M = MeV2_to_Tesla * abs(compute_magnetization(T, m, g, e_val, B, Omega))
        M_vals.append(M)
        #print(f"T = {T:.2f} MeV, M = {M:.6e}")

    # Plot the magnetization vs. temperature.
    plt.figure(figsize=(8,6))
    plt.plot(T_vals, M_vals, marker='o', linestyle='-', label='Magnetization')
    plt.grid(True)
    plt.xlabel("Temperature (MeV)", fontsize=14)
    plt.ylabel("Magnetic Field (T)", fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(1e0, 1e18)
    plt.gca().invert_xaxis()
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

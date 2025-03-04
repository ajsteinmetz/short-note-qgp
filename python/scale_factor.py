import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Constants
T0 = 2.35e-4         # Today's CMB temperature in eV
B0_low = 1e-20       # Lower end of today's intergalactic magnetic field (Tesla)
B0_high = 1e-12      # Upper end (Tesla)
B_critical = 4.41e9  # Schwinger critical field in Tesla

# Quark magnetization constant at 300 MeV
Bq_300 = 6.958e15    # Tesla at T = 300 MeV

# Temperature range in MeV; we'll use 500 points between 50 and 500 MeV.
T_MeV = np.linspace(50, 500, 1000)
# Convert temperature from MeV to eV (1 MeV = 1e6 eV)
T_eV = T_MeV * 1e6

# Compute the primordial magnetic field using the scaling law:
# B_PMF(T) = B0 * (T^2)/(T0^2)
B_low = B0_low * (T_eV**2) / (T0**2)
B_high = B0_high * (T_eV**2) / (T0**2)

# Split the temperature range for the quark magnetization:
# 1. For T >= 170 MeV: use the standard B_q formula.
mask_orig = T_MeV >= 170
T_quark_orig = T_MeV[mask_orig]
B_q_orig = Bq_300 / (300**3) * T_quark_orig**3

# 2. For the transition region 150 MeV <= T < 170 MeV:
mask_trans = (T_MeV >= 150) & (T_MeV < 170)
T_quark_trans = T_MeV[mask_trans]
# f(T) = (T - 150)/20, which linearly scales from 0 at 150 MeV to 1 at 170 MeV.
f = (T_quark_trans - 150) / 20
B_q_trans = (Bq_300 / (300**3) * T_quark_trans**3) * f

# Compute the electron magnetization over the full temperature range:
# B_e(T) = 2.15 * B_q(T) where B_q is computed for the full range using the original formula.
B_e = 2.15 * (Bq_300 / (300**3) * T_MeV**3)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Shade the region between the lower and upper bounds for the PMF
ax.fill_between(T_MeV, B_low, B_high, color='lightgrey', alpha=1.0, label='Allowed PMF range')

# Plot the PMF curves
ax.plot(T_MeV, B_low, label=f'$B_0 = {B0_low}$ T', color='black')
ax.plot(T_MeV, B_high, label=f'$B_0 = {B0_high}$ T', color='black')

# Plot the quark magnetization curves:
# - The original curve for T >= 170 MeV (red solid line)
ax.plot(T_quark_orig, B_q_orig, label=r'$B_q(T)$ for $T \geq 170$ MeV', 
        color='red', linewidth=2)
# - The transition curve for 150 MeV <= T < 170 MeV (red solid line)
ax.plot(T_quark_trans, B_q_trans, label=r'$B_q(T)\times\frac{T-150}{20}$ for $150\leq T < 170$ MeV', 
        color='red', linestyle='--', linewidth=2)

# Plot the electron magnetization curve (for all T)
ax.plot(T_MeV, B_e, label=r'$B_e(T) = 2.15\,B_q(T)$', color='blue', linewidth=2)

# Plot the horizontal line for the Schwinger critical field
ax.axhline(y=B_critical, color='black', linestyle='-.', linewidth=1.5,
           label=f'$B_C = {B_critical:.2e}$ T')

# Add a vertical dashed red line at T = 150 MeV
ax.axvline(x=150, color='green', linestyle='-', linewidth=1.5, label=r'$T = 150\,\mathrm{MeV}$')

# Set labels
ax.set_xlabel("Temperature (MeV)", fontsize=14)
ax.set_ylabel("Magnetic Field (T)", fontsize=14)

# Set y-axis to log scale and invert the x-axis so that the universe gets cooler from left to right
ax.set_yscale('log')
ax.invert_xaxis()
ax.set_xlim(500, 50)

# Add gridlines
ax.grid(which='major', linestyle='--', linewidth=0.75, alpha=1.0)
ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=100))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

plt.show()

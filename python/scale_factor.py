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
T_MeV = np.linspace(50, 500, 500)
# Convert temperature from MeV to eV (1 MeV = 1e6 eV)
T_eV = T_MeV * 1e6

# Compute the primordial magnetic field using the scaling law:
# B_PMF(T) = B0 * (T^2)/(T0^2)
B_low = B0_low * (T_eV**2) / (T0**2)
B_high = B0_high * (T_eV**2) / (T0**2)

# Compute the quark magnetization as a function of temperature for T > 150 MeV:
# B_q(T) = Bq_300/(300^3) * T^3
mask = T_MeV > 150
T_quark = T_MeV[mask]
B_q = Bq_300 / (300**3) * T_quark**3

# Compute the electron magnetization over the full temperature range:
# B_e(T) = 2.15 * B_q(T) where B_q is computed for the full range
B_e = 2.15 * (Bq_300 / (300**3) * T_MeV**3)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Shade the region between the lower and upper bounds for the PMF, quark magnetization, and electron magnetization
#ax.fill_between(T_quark, 0, B_q, color='red', alpha=0.10, label='Quark magnetization range')
ax.fill_between(T_MeV, B_low, B_high, color='lightgrey', alpha=1.0, label='Allowed PMF range')
# ax.fill_between(T_MeV, 0, B_e, color='blue', alpha=0.05, label='Electron magnetization range')

# Plot the PMF curves
ax.plot(T_MeV, B_low, label=f'$B_0 = {B0_low}$ T', color='black')
ax.plot(T_MeV, B_high, label=f'$B_0 = {B0_high}$ T', color='black')

# Plot the quark magnetization curve (only for T > 150 MeV)
ax.plot(T_quark, B_q, label=r'$B_q(T) = B_q(300\,\mathrm{MeV}) \left(\frac{T}{300\,\mathrm{MeV}}\right)^3$', 
        color='red', linewidth=2)

# Plot the electron magnetization curve (for all T)
ax.plot(T_MeV, B_e, label=r'$B_e(T) = 2.15\,B_q(T)$', color='blue', linewidth=2)

# Plot the horizontal line for the Schwinger critical field
ax.axhline(y=B_critical, color='black', linestyle='-.', linewidth=1.5,
           label=f'$B_C = {B_critical:.2e}$ T')

# Add a vertical dashed black line at T = 150 MeV
ax.axvline(x=150, color='red', linestyle='--', linewidth=1.5, label=r'$T = 150\,\mathrm{MeV}$')

# Set labels and title
ax.set_xlabel("Temperature (MeV)", fontsize=14)
ax.set_ylabel("Magnetic Field (T)", fontsize=14)
# ax.set_title("Magnetic Fields vs. Temperature", fontsize=16)

# Set y-axis to log scale and invert the x-axis so that the universe gets cooler from left to right
ax.set_yscale('log')
ax.invert_xaxis()
ax.set_xlim(500, 50)

# Add major gridlines
ax.grid(which='major', linestyle='--', linewidth=0.75, alpha=1.0)

# Explicitly set dashed minor ticks and gridlines for the y-axis using LogLocator
ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=100))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
#ax.grid(which='minor', axis='y', linestyle='--', linewidth=0.5, alpha=1.0)

# Set minor ticks for the x-axis and add minor gridlines for x-axis
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
#ax.grid(which='minor', axis='x', linestyle='--', linewidth=0.5, alpha=1.0)

# Commented out the legend
# ax.legend(fontsize=12)

plt.show()

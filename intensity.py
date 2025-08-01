
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
import mendeleev
import re
from basefunctions import *

# Constants
Na = const.Avogadro
r_e = (const.e**2)/(4*np.pi*const.epsilon_0*const.m_e*const.c**2)

# Parameters
n1=1.05
n2 = 1.00
N_photons = 1000
pore_fraction = 0.985
pore_radius = 1e-4 #1e-4 to 1e-7
foam_height = 1700e-4 #cm
foam_diameter = 900e-4 #cm
foam_radius = foam_diameter/2
foam_volume = np.pi*((foam_diameter)/2)**2*foam_height
pore_volume = 4/3*(np.pi)*pore_radius**3
N_pores = int(pore_fraction*(foam_volume)/pore_volume)

def parse_compound(formula):
    """Parse chemical formula into elements and counts"""
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = []
    counts = []
    for elem, count in matches:
        elements.append(elem)
        counts.append(int(count) if count else 1)
    return elements, counts

#returns only WEIGHTED-MU
def calculate_mu_compound(formula, rho, energy_keV):
    """Calculate attenuation coefficient for a compound"""
    elements, counts = parse_compound(formula)
    energy_eV = energy_keV * 1000  # Convert to eV
    
    # Initialize arrays
    mu_total = np.zeros_like(energy_eV)
    total_mass = 0
    
    # Calculate molar mass and prepare weights
    molar_masses = []
    for elem in elements:
        molar_masses.append(mendeleev.element(elem).mass)
    molar_masses = np.array(molar_masses)
    counts = np.array(counts)
    total_mass = np.sum(molar_masses * counts)
    
    # Weighted sum of mu/rho
    weighted_mu_rho = 0
    for elem, count, mm in zip(elements, counts, molar_masses):
        energy, f1, f2 = load_file(elem)
        f2_interp = np.interp(energy_eV, energy, f2)
        lambda1 = energy_to_wavelength(energy_keV)
        mu_rho_elem = (Na * 2 * r_e * lambda1 * f2_interp)/mm
        weighted_mu_rho += (count * mm/total_mass) * mu_rho_elem
    
    return weighted_mu_rho * rho


def simulate_horizontal_displacement(N_photons, foam_radius, pore_radius, pore_fraction):
    avg_pores = avg_interactions(pore_fraction, foam_height, pore_radius)
    displacements = np.zeros(N_photons)
    for i in range(N_photons):
        dx_tot = 0.0
        theta_i = 0  # Initial angle (normal incidence)
        for _ in range(avg_pores):
            d = np.random.uniform(-pore_radius, pore_radius)
            theta_i = angle_of_incidence(d, pore_radius)
            theta_t = snells_law(theta_i, n1, n2)
            dx = lateral_displacement(theta_i, theta_t, pore_radius)
            dx_tot += dx
            theta_i = 2*theta_t - theta_i  # Update for next interaction
        displacements[i] = dx_tot
    return displacements

def calculate_transmission(energy_keV, formula, rho, plot_hist=True):
    # Get horizontal displacements
    delta_x = simulate_horizontal_displacement(N_photons, foam_radius, pore_radius, pore_fraction)
    
    # Calculate path lengths
    path_lengths = np.sqrt(foam_height**2 + delta_x**2)
    
    # Get attenuation coefficient
    mu = calculate_mu_compound(formula, rho, energy_keV)

    # Calculate transmissions
    transmissions = np.exp(-mu * path_lengths)
    
    #calculate percent differences
    pure_transmission = np.exp(-1*mu*foam_height)
    mean_scattering_transmission = np.mean(transmissions)
    print(f'Pure attenuation would result in an intensity of {pure_transmission}')
    print(f'The average intensity consider scattering yield {mean_scattering_transmission}')
    print(f'The percent difference caused by scattering is {(np.abs(pure_transmission-mean_scattering_transmission)/pure_transmission)*100}')
    # Plot histogram
    if plot_hist:
        plt.hist(transmissions, 
             bins='auto',
             color='royalblue',
             edgecolor='navy',
             alpha=0.8)
        from matplotlib.ticker import ScalarFormatter
        # Force proper decimal formatting
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=False))
        ax.ticklabel_format(style='plain', axis='both')  # Disable scientific notation
        
        # Formatting
        plt.xlabel('Transmission Value', fontsize=12)
        plt.ylabel('Number of Photons', fontsize=12)
        plt.title('Photon Transmission Distribution with Scattering', fontsize=14)
        plt.grid(True, linestyle=':', alpha=0.4)
        
        plt.tight_layout()
        plt.show()

energy_keV = 2.0  # keV
formula = "Si5O10Ti"   # Compound formula
rho = 2.33       # g/cm³
mu = calculate_mu_compound(formula, rho, energy_keV)
transmissions = calculate_transmission(energy_keV, formula, rho)

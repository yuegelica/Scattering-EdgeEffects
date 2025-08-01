import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
import mendeleev
import re
from basefunctions import *

# Constants
Na = const.Avogadro
r_e = (const.e**2)/(4*np.pi*const.epsilon_0*const.m_e*const.c**2)

def load_file(filename):
    data = np.loadtxt(filename, skiprows=1)
    energy = data[:,0]  # in eV
    f1 = data[:,1]
    f2 = data[:,2]
    lambda1 = energy_to_wavelength(energy)  # in meters
    return energy, f1, f2, lambda1

def parse_compounds(formula):
    format = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(format, formula)
    elements = []
    counts = []
    for (elem, count) in matches:
        elements.append(elem)
        if count:
            counts.append(int(count))
        else:
            counts.append(1)
    return elements, counts


def weighted_ma(formula, rho, x):
    elements, counts = parse_compounds(formula)
    
    mu_rho_list = []
    energy_common = None
    
    for i, elem in enumerate(elements):
        filename = f"{elem.lower()}.txt"
        energy, f1, f2, lambda1 = load_file(filename)
        rho_si = 1000.0  # kg/m^3 for mass attenuation calculation
        mu = attenuation(f2, rho_si, lambda1, elem)
        mu_rho = mu / rho_si  # m^2/kg
        
        if energy_common is None:
            energy_common = energy # use first element's energy grid
            mu_rho_list.append(mu_rho)
        else:
            # interpolate onto common energy grid
            mu_rho_interp = np.interp(energy_common, energy, mu_rho)
            mu_rho_list.append(mu_rho_interp)
    
    # Stack the mu/rho arrays
    mu_rho_array = np.vstack(mu_rho_list)
    
    # Calculate molar masses and weight fractions
    molar_masses = np.array([mendeleev.element(e).mass for e in elements])
    counts = np.array(counts)
    total_mass = np.sum(molar_masses * counts)
    weight_fractions = (molar_masses * counts) / total_mass
    
    # Calculate weighted average of mu/rho
    weighted_mu_rho = np.zeros_like(mu_rho_array[0])
    for i in range(len(weight_fractions)):
        weighted_mu_rho += weight_fractions[i] * mu_rho_array[i]
    
    rho_si = rho * 1000.0  # kg/m^3
    
    mu = weighted_mu_rho * rho_si  # m^-1
    
    x_si = x * 0.01  # meters
    
    transmission = calculate_transmission(mu, x_si)
    
    return energy_common, weighted_mu_rho, transmission

def plot_transmission(photon_energy, transmission, reference_data=None):
    """Plot transmission vs photon energy"""
    plt.figure(figsize=(12, 8))
    plt.xlabel("Photon Energy (eV)")
    plt.ylabel("Transmission")
    plt.title("X-Ray Transmission of Si5O10Ti")
    plt.xlim(1000, 10000)
    plt.ylim(0.92, 1)
    
    # Plot calculated transmission
    plt.plot(photon_energy, transmission, 'b-', linewidth=2, label='Calculated')
    
    # Plot reference data if provided
    if reference_data is not None:
        ref_energy, ref_transmission = reference_data
        plt.plot(ref_energy, ref_transmission, 'r--', linewidth=2, label='Reference')
    
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

def calculate_and_plot_transmission():
    # Parameters
    formula = "Si5O10Ti"
    rho = 2.33  # g/cm^3 (density of material)
    x = 1e-5    # cm (0.1 microns thickness)
    
    # Reference data
    reference_energy = np.array([
        1000.0, 1090.0, 1180.0, 1270.0, 1360.0, 1450.0, 1540.0, 1630.0, 1720.0, 1810.0,
        1900.0, 1990.0, 2080.0, 2170.0, 2260.0, 2350.0, 2440.0, 2530.0, 2620.0, 2710.0,
        2800.0, 2890.0, 2980.0, 3070.0, 3160.0, 3250.0, 3340.0, 3430.0, 3520.0, 3610.0,
        3700.0, 3790.0, 3880.0, 3970.0, 4060.0, 4150.0, 4240.0, 4330.0, 4420.0, 4510.0,
        4600.0, 4690.0, 4780.0, 4870.0, 4960.0, 5050.0, 5140.0, 5230.0, 5320.0, 5410.0,
        5500.0, 5590.0, 5680.0, 5770.0, 5860.0, 5950.0, 6040.0, 6130.0, 6220.0, 6310.0,
        6400.0, 6490.0, 6580.0, 6670.0, 6760.0, 6850.0, 6940.0, 7030.0, 7120.0, 7210.0,
        7300.0, 7390.0, 7480.0, 7570.0, 7660.0, 7750.0, 7840.0, 7930.0, 8020.0, 8110.0,
        8200.0, 8290.0, 8380.0, 8470.0, 8560.0, 8650.0, 8740.0, 8830.0, 8920.0, 9010.0,
        9100.0, 9190.0, 9280.0, 9370.0, 9460.0, 9550.0, 9640.0, 9730.0, 9820.0, 9910.0,
        10000.0
    ])

    reference_transmission = np.array([
        0.92036, 0.93612, 0.94781, 0.95678, 0.96385, 0.96944, 0.97397, 0.97767, 0.98069, 0.98321,
        0.95919, 0.96357, 0.96733, 0.97064, 0.97353, 0.97606, 0.97829, 0.98025, 0.98199, 0.98353,
        0.98491, 0.98614, 0.98724, 0.98823, 0.98912, 0.98992, 0.99065, 0.99131, 0.99191, 0.99245,
        0.99295, 0.99341, 0.99383, 0.99421, 0.99457, 0.99489, 0.99519, 0.99547, 0.99572, 0.99596,
        0.99618, 0.99639, 0.99658, 0.99675, 0.99692, 0.99514, 0.99536, 0.99558, 0.99577, 0.99596,
        0.99614, 0.99630, 0.99646, 0.99661, 0.99675, 0.99688, 0.99701, 0.99713, 0.99724, 0.99735,
        0.99745, 0.99755, 0.99764, 0.99773, 0.99781, 0.99789, 0.99797, 0.99804, 0.99811, 0.99817,
        0.99823, 0.99829, 0.99835, 0.99841, 0.99846, 0.99851, 0.99856, 0.99860, 0.99865, 0.99869,
        0.99873, 0.99877, 0.99880, 0.99884, 0.99887, 0.99891, 0.99894, 0.99897, 0.99900, 0.99903,
        0.99905, 0.99908, 0.99911, 0.99913, 0.99915, 0.99918, 0.99920, 0.99922, 0.99924, 0.99926,
        0.99928
    ])

    # Calculate transmission using the compound formula
    energy_common, weighted_mu_rho, transmission = weighted_ma(formula, rho, x)
    
    # Plot the result with reference data
    plot_transmission(energy_common, transmission, (reference_energy, reference_transmission))
    
    return energy_common, transmission

if __name__ == "__main__":
    energy, transmission = calculate_and_plot_transmission()

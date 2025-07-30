import numpy as np
import mendeleev
import scipy.constants as const
import re

### Constants
Na = const.Avogadro
h = const.h
hbar = const.hbar
c = const.c
hc = h * c
r_e = (const.e**2) / (4 * const.pi * const.epsilon_0 * const.m_e * const.c**2)

### Conversion functions
def energy_to_wavelength(energy):
    energy_joules = energy * const.e  # Convert eV to J
    return hc / energy_joules
#wavelength will be in meters for these calculations
def wavelength_to_energy(wavelength):
    energy_joules =  hc / wavelength
    return energy_joules / const.e

#angular frequency
def energy_to_frequency(energy):
    energy_joules = energy * const.e
    return energy_joules / h

def frequency_to_energy(frequency):
    return h*frequency / const.e

def frequency_to_omega(frequency):
    return 2 * np.pi * frequency

def wavelength_to_wavenumber(wavelength):
    return 2 * np.pi / wavelength

def wavenumber_to_momentum(k):
    return hbar *k

###############################################################################
# Load data from Henke file, skip the first row as those are the labels (str)
def load_file(filename: str):
    energy = np.loadtxt(filename, usecols=0, skiprows=1)# in eV
    f1 = np.loadtxt(filename, usecols=1, skiprows=1)
    f2 = np.loadtxt(filename, usecols=2, skiprows=1)
    lambda1 = energy_to_wavelength(energy)# in meters
    return energy, f1, f2, lambda1

# Attenuation coefficient mu
def attenuation(f2, rho, lambda1, element):
    elem = mendeleev.element(element)
    Ma = elem.mass  # g/mol
    mu = (rho * Na * 2 * r_e * lambda1 * f2) / Ma
    return mu

# Mass attenuation coefficient mu/rho
def mass_attenuation(mu, rho):
    mass_attenuation = mu/rho
    return mass_attenuation

# Transmission, account for absorbtion only (f2)
def calculate_transmission(mu, x):
    return np.exp(-mu * x)

# Find delta (the real refractive index) from f1
def delta(f1, rho, lambda1):
    delta = (rho*lambda1**2)/(2*np.pi*Na*f1)
    return delta

#Find the index of refraction n2 from delta
def n2(delta):
    n2 = 1- delta
    return n2

# Parses compounds (elements and counts) to use the correct files and values
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

#for graphing transmission (transmissionplots.py)
def weighted_ma_graph(formula, rho, x):
    elements, counts = parse_compounds(formula)

    mu_rho_list = []
    energy_common = None

    for elem in elements:
        filename = f"{elem.lower()}.txt"
        energy, f1, f2, lambda1 = load_file(filename)

        mu = attenuation(f2, f1, lambda1, elem)
        mu_rho = mass_attenuation(mu, 1.0)

        if energy_common is None:
            energy_common = energy  # use first element's energy grid but should try to make a uniform one in the future
            mu_rho_list.append(mu_rho)
        else:
            # interpolate onto common_energy grid
            mu_rho_interp = np.interp(energy_common, energy, mu_rho)
            mu_rho_list.append(mu_rho_interp)

    mu_rho_array = np.vstack(mu_rho_list)

    molar_masses = np.array([mendeleev.element(e).mass for e in elements])
    counts = np.array(counts)
    total_mass = np.sum(molar_masses * counts)
    weight_fractions = (molar_masses * counts) / total_mass

    weighted_mu_rho = np.zeros_like(mu_rho_array[0])
    for i in range(len(weight_fractions)):
        weighted_mu_rho += weight_fractions[i] * mu_rho_array[i]

    mu = weighted_mu_rho * rho
    transmission = np.exp(-mu * x)

    return energy_common, weighted_mu_rho

# The incident angle theta_i depending on the distance from the center it hits the pore
def angle_of_incidence(d, r):
    return np.arcsin(d/r)

# Find theta_refracted
def snells_law(theta_i, n1, n2):
    sin_theta_i = np.sin(theta_i)
    sin_theta_t = (n1/n2) * sin_theta_i

    # Clip to valid range [-1, 1] to prevent arcsin errors
    sin_theta_t_clipped = np.clip(sin_theta_t, -0.999999, 0.999999)
    
    if np.abs(sin_theta_t_clipped) > 1:
        return theta_i  # reflect, i.e., angle stays the same
    else:
        return np.arcsin(sin_theta_t_clipped)

#Find how much a photon gets horizontally displaced ore after interacting with and being refracted through each pore
def lateral_displacement(theta_i, theta_t, r):
    return (2*(r**2)-2*(r**2)*np.cos(np.pi-2*theta_t))**0.5*np.sin(theta_t-theta_i)
    # return r * 2 *  ( np.tan(theta_t) - np.tan(theta_i))

# Find the average pores a photon will hit from the mean free path
def calc_avg_interactions(foam_length, pore_radius, pore_fraction):
    mean_free_path = (4/3)*pore_radius*(1-pore_fraction)/pore_fraction
    return int(foam_length/mean_free_path)


###############################################################
# Useful functions but not used directly in this sumulation #
def critical_angle(f1, rho,lambda1, element):
    elem = mendeleev.element(element)
    Ma = elem.mass
    theta_c = np.sqrt((rho*Ma*r_e*lambda1**2)/np.pi(Na))*f1
    return theta_c
# Optical_depth is absorption length
def optical_depth(mu, rho):
    optical_depth = rho/mu
    return optical_depth
#Porosity given density of solid and aerogel
def calc_porosity(rho_solid, rho_eff):
    return 1- (rho_eff/rho_solid)

# Fresnel Equations to calculate how much specular reflection you expect 
def complex_reflection(n1, n2, theta_i, theta_t):
    theta_t = np.arcsin(n1*np.sin(theta_i/n2))
    r_s = np.abs((n1*np.cos(theta_i)-n2*np.cos(theta_t))/(n1*np.cos(theta_i)+n2*np.cos(theta_t)))**2
    r_p = np.abs((n2*np.cos(theta_i)-n1*np.cos(theta_t))/(n2*np.cos(theta_i)+n1*np.cos(theta_t)))**2
    t_s = 1-r_s
    t_p = 1-r_p
    return r_s, r_p, t_s, t_p

# Unpolarized natural right
def r_effective(r_s, r_p):
    r_eff = 0.5*(r_s+r_p)
    return r_eff

#Returns the distance a photons gets refracted through a rectangular slab
def refraction(theta_i,foam_height):
    d = foam_height/ (np.tan((90-theta_i)))
    return d 

def snells_law(theta_i, n1, n2):
    sin_theta_i = np.sin(theta_i)
    sin_theta_t = (n1 / n2) * sin_theta_i
#special case where theta_i = theta_t = 0, so find specular reflection to show it is negligible if x-rays go thorugh a foam
def normal_incidence(n1, n2):
    reflectance = np.abs((n1-n2)/(n1+n2))**2
    return reflectance

# Find beta, (imaginary component of refractive index) from f2
def beta(f2, lambda1, element, rho):
    elem= mendeleev.element(element)
    Ma = elem.mass
    beta = (rho*Ma*r_e*f2*lambda1**2)/(2*np.pi*Na)
    return beta

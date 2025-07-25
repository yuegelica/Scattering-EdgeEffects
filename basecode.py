import numpy as np
import mendeleev
import scipy.constants as const
import re

# Constants
Na = const.Avogadro
h = const.h
hbar = const.hbar
c = const.c
hc = h * c
r_e = (const.e**2) / (4 * const.pi * const.epsilon_0 * const.m_e * const.c**2)

# Conversion functions
def energy_to_wavelength(energy):
    energy_joules = energy * const.e  #Convert eV to J
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

# Load data from Henke file
def load_file(filename: str):
    energy = np.loadtxt(filename, usecols=0, skiprows=1)# in eV
    f1 = np.loadtxt(filename, usecols=1, skiprows=1)
    f2 = np.loadtxt(filename, usecols=2, skiprows=1)
    lambda1 = energy_to_wavelength(energy)# in m
    return energy, f1, f2, lambda1

# Attenuation coefficient mu
def attenuation(f2, rho, lambda1, element):
    elem = mendeleev.element(element)
    Ma = elem.mass  # g/mol
    mu = (rho * Na * 2 * r_e * lambda1 * f2) / Ma
    return mu

#Mass attenuation coefficient mu/rho
def mass_attenuation(mu, rho):
    mass_attenuation = mu/rho
    return mass_attenuation

# Transmission, account for absorbtion only (f2)
def calculate_transmission(mu, x):
    return np.exp(-mu * x)

#optical_depth is absorbtion length
def optical_depth(mu, rho):
    optical_depth = rho/mu
    return optical_depth

#Fresnel Equations to calculate how much specular reflection you expect 
def complex_reflection(n1, n2, theta_i, theta_t):
    theta_t = np.arcsin(n1*np.sin(theta_i/n2))
    r_s = np.abs((n1*np.cos(theta_i)-n2*np.cos(theta_t))/(n1*np.cos(theta_i)+n2*np.cos(theta_t)))**2
    r_p = np.abs((n2*np.cos(theta_i)-n1*np.cos(theta_t))/(n2*np.cos(theta_i)+n1*np.cos(theta_t)))**2
    t_s = 1-r_s
    t_p = 1-r_p
    return r_s, r_p, t_s, t_p

#unpolarized natural right
def r_effective(r_s, r_p):
    r_eff = 0.5*(r_s+r_p)
    return r_eff

#find delta
def delta(f1, rho, lambda1):
    delta = (rho*lambda1**2)/(2*np.pi*Na*f1)
    return delta
#find n2 from delta
def n2(delta):
    n2 = 1- delta
    return n2

#special case where theta_i = theta_t = 0, so find specular reflection to show it is negligible if x-rays go thorugh a foam
def normal_incidence(n1, n2):
    reflectance = np.abs((n1-n2)/(n1+n2))**2
    return reflectance

# beta = imagenary refractive index
def beta(f2, lambda1, element, rho):
    elem= mendeleev.element(element)
    Ma = elem.mass
    beta = (rho*Ma*r_e*f2*lambda1**2)/(2*np.pi*Na)
    return beta

#use f1 to find othre quantities
#actual scatter factor, plug into fresnel's 
def critical_angle(f1, rho,lambda1, element):
    elem = mendeleev.element(element)
    Ma = elem.mass
    theta_c = np.sqrt((rho*Ma*r_e*lambda1**2)/np.pi(Na))*f1
    return theta_c


#returns amount deflected horizontally through rectangular foam
def refraction(theta_i,foam_length):
    d = foam_length/ (np.tan((90-theta_i)))
    return d 

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

    return energy_common, transmission

#Now, consider statistical model, lots of photons come in starting parallel, given porosity of foam and nmber of pores, we calculate the distribution of horizontal refractions. 
#Dictated by graph, and probability of location photon is hitting is uniform, the particular angle of incidence probability follows the curve
#look in numpy random for sample

def angle_of_incidence(d, r):
    return np.arcsin(d/r)

def snells_law(theta_i, n1, n2):
    # Calculate sin(theta_t) with clipping
    sin_theta_i = np.sin(theta_i)
    sin_theta_t = (n1/n2) * sin_theta_i

    sin_theta_t_clipped = np.clip(sin_theta_t, -0.999999, 0.999999)
    
    return np.arcsin(sin_theta_t_clipped)

def lateral_displacement(theta_i, theta_t, r):
    return (2*(r**2)-2*(r**2)*np.cos(np.pi-2*theta_t))**0.5*np.sin(theta_t-theta_i)

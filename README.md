# Scattering-EdgeEffects
## Motivation
The project aims to investigate whether the inconsistent density measurements seen at the edges of an aerogel foam sample are a consequence of some sort of scattering mechanism (whether due to the porosity, surface mechanisms like refraction, or something else entirely). 
## Table of Contents
1. [Motivation](#motivation)  
2. [Instructions to Run](#instructions-to-run)  
   - [Required Installations](#required-installations)  
   - [Cloning the Repository](#cloning-the-repository)  
   - [Dependencies](#dependencies)  
3. [Running the Simulation](#running-the-simulation)  
   - [CCD Distribution](#ccd-distribution)  
   - [Transmission Plot](#transmission-plot)  
   - [Horizontal Displacement](#horizontal-displacement)  
   - [Intensity Distribution](#intensity-distribution)  
4. [Usage](#usage)  
5. [Customizations](#customizations)  
6. [Recreating Simulation Results](#recreating-simulation-results)  
7. [Paper and Presentation](#paper-and-presentation)  
8. [Contact](#contact)  
9. [Useful Resources](#useful-resources)
## Instructions to Run
### Required Installations
Installing a virtual envionment is recommended:
```bash
python3 -m venv myenv
```
```bash
source myenv/bin/activate
```
Install the following Python libraries:

```bash
pip install numpy scipy matplotlib mendeleev
```
### Cloning the Repository
```bash
git clone git@github.com:yuegelica/Scattering-EdgeEffects.git
cd Scattering-EdgeEffects
```
### Dependencies
Install the required packages listed in requirements.txt:
```bash
pip install -r requirements.txt
```

## Running the Simulation
There are 4 different plots this project can produce. The `basecode.py` file contains all the conversions and functions needed for the simulations to run.
### CCD Distribution
The main plot can be generated in the `main.py` file. It plots the combined intensity values at each CCD pixel, the zoomed in left and right edges of the foam, as well as the values of the difference in transmission.
### Transmission Plot
To generate a Transmimssion plot identical to ones found on https://henke.lbl.gov/, run the `transmissionplot.py` file. You will need the formula of the compound, path length, and density based off the composition of the material. 
### Horizontal Displacement
To generate a histogram of the horizontal displacement distribution, run the `horizontaldisplacement.py` file. 
### Intensity Distribution
To generate a histogram of the intensity distribution, run the `intensity.py` file. 
*Filenames to be changed and files to be combined*
[Back to Top](#table-of-contents)
## Customizations
You are free to modify
1. Formula of Compound/Element
2. porosity (pore_fraction)
3. Refractive index of the Foam (n2)
4. Pore Radius (pore_radius)
5. Number of photons (N_photons)
## Recreating Simulation Results
Refer to `recreate_plots.py`
## Paper and Presentation
1. To build paper, change directories into `paper/`
```
cd paper/
```
2. To build presentation, change directories into `presentation/`
```
cd presentation/
```
## Contact
If you have any questions or suggestions, feel free to contact me at yuegelicay@gmail.com.

### Useful Resources


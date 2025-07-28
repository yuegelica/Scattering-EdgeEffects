# Scattering-EdgeEffects
The project aims to investigate whether the inconsistent measurements seen at the edges of an aerogel foam sample are a consequence of some sort of scattering mechanism (whether due to the porosity, surface mechanisms like refraction, or something else entirely). 
## Table of Contents

## Motivation

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
### Dependencies
Install the required packages listed in requirements.txt:
```bash
pip install -r requirements.txt
```
### Cloning the Repository
```bash
git clone git@github.com:yuegelica/Scattering-EdgeEffects.git
```

## Running the Simulation
There are 4 different plots this project can produce. The basecode.py file contains all the conversions and functions needed for the simulations to fun.
### CCD Distribution
The main plot can be generated in the main.py file. It plots the combined intensity values at each CCD pixel.There are three subplots that will be generated.
### Transmission Plot
To generate a Transmimssion plot identical to ones found on https://henke.lbl.gov/, run the transmissionplot.py file.
### Horizontal Displacement
To generate a histogram of the horizontal displacement distribution, run the horizontaldisplacement.py file.
### Intensity Distribution
To generate a histogram of the intensity distribution, run the intensity.py file.
## Usage

### Example Run

## Recreating Simulation Results
Refer to `recreate_plots.py`
## Paper and Presentation

## Contact
If you have any questions or suggestions, feel free to contact me at yuegelicay@gmail.com.

### Useful Resources


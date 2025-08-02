# Scattering-EdgeEffects
## Motivation
The project aims to investigate whether the inconsistent density measurements seen at the edges of an aerogel foam sample are a consequence of some sort of scattering mechanism (whether due to the porosity, surface mechanisms like refraction, or something else entirely). 
## Table of Contents
1. [Motivation](#motivation)
2. [Setup and Installation](#setup-and-installation)
   - [Required Installations](#required-installations)
   - [Cloning the Repository](#cloning-the-repository)
   - [Dependencies](#dependencies)
3. [Running Simulations](#running-simulations)
   - [Main CCD Distribution Plots](#1-generate-main-ccd-distribution-plots)
   - [Comprehensive Summary Plots](#2-generate-comprehensive-summary-plots)
   - [X-ray Transmission Plots](#3-generate-x-ray-transmission-plots)
   - [Individual Analysis Plots](#4-generate-individual-analysis-plots)
   - [Relationship Summary Plots](#5-generate-relationship-summary-plots)
   - [Recreate All Results](#6-recreate-all-simulation-results)
4. [Customizations](#customizations)
5. [File Structure and Workflow](#file-structure-and-workflow)
   - [Core Simulation Files](#core-simulation-files)
   - [Expected Outputs](#expected-outputs)
6. [Paper and Presentation](#paper-and-presentation)
   - [Building the Paper](#building-the-paper)
   - [Building the Presentation](#building-the-presentation)
7. [Contact](#contact)
8. [Useful Resources](#useful-resources)
## Setup and Installation
### Required Installations
Installing a virtual envionment is recommended:
```bash
python3 -m venv myenv
```
```bash
source myenv/bin/activate
```
# Cloning the Repository:
```bash
git clone git@github.com:yuegelica/Scattering-EdgeEffects.git
cd Scattering-EdgeEffects
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
## Running Simulations
There are multiple different plots this project can produce. The `basecode.py` file contains all the conversions and functions needed for the simulations to run.

### 1. Generate main CCD distribution plots

Run the main simulation that shows combined intensity values at each CCD pixel:
```
python main.py
```

This creates plots showing:
- Combined intensity values at each CCD pixel: CCD_Overlay_new.png
- Zoomed view of left edge effects: zoomed_left_edge_new.png
- Zoomed view of right edge effects: zoomed_right_edge_new.png
- Transmission difference values: transmission_diff.png
### 2. Generate comprehensive summary plots

Create six summary subplots:
```
python mainraw.py
```

This produces:
- Input beam profile: inputbeamprofile.png
- CCD sensor with pure attenuation: pureattenuationccd.png
- Horizontal displacement distribution: horizontaldisplacements.png
- Transmission distribution comparison: scattering_transmission_comparison.png
- CCD sensor with scattering effects: scatteredccd.png
- Final position histogram: final_ray_pos.png

### 3. Generate X-ray transmission plots

Create transmission vs energy plots (similar to Henke database): 
```
python transmissionplot.py
```

This generates transmission curves for your material composition. For "Si5O10Ti", you would get the file energy_vs_transmission.png

### 4. Generate individual analysis plots

Create horizontal displacement histogram:
```
python horizontaldisplacement.py
```
And this would create a file called lateral_displacement.png
Create intensity distribution histogram:
```
python intensity.py
```

### 5. Generate relationship summary plots

Create trend analysis plots showing effects of pore radius and porosity:
```
python summaryplots.py
```

This produces four plots summarizing:
- Pore radius vs transmission: dxvsporeradius.png
- Pore radius vs horizontal displacement: dxvsporositynew.png
- Porosity vs transmission: transmissionvsporosity1.png
- Porosity vs horizontal displacement: transmissonvspore1.png

### 6. Recreate all simulation results

To use the default variables and reproduce the results:
```
python recreate_plots.py
```
[Back to Top](#table-of-contents)
## Customizations
You can modify these parameters in the relevant scripts:
1. **Material composition**: Formula of compound/element (e.g., "Si5O10Ti")
2. **Porosity**: `pore_fraction` parameter (0.0 to 1.0)
3. **Refractive index**: `n2` parameter for foam material
4. **Pore size**: `pore_radius` parameter in meters
5. **Simulation size**: `N_photons` parameter for statistical accuracy

## File structure and workflow

The typical workflow is:

1. Edit simulation parameters in the Python files
2. Run the simulations to generate data and plots:
   ```
   python main.py
   python mainraw.py  
   python transmissionplot.py
   python summaryplots.py
   ```
3. Update the paper with new results:
   ```
   cd paper/
   pdflatex main.tex
   bibtex main
   pdflatex main.tex
   pdflatex main.tex
   ```
4. Update the presentation:
   ```
   cd presentation/
   pdflatex CleanEasy.tex
   ```

## Core simulation files

- `basecode.py` - Contains all conversion functions and core simulation steps
- `main.py` - Main CCD distribution analysis
- `mainraw.py` - Comprehensive six-plot summary
- `transmissionplot.py` - X-ray transmission calculations
- `horizontaldisplacement.py` - Displacement analysis
- `intensity.py` - Intensity distribution analysis
- `summaryplots.py` - Parameter relationship plots
- `recreate_plots.py` - Reproduce all results with consistent parameters

## Expected outputs

After running all simulations, you should have:
- Multiple PNG plots showing various aspects and steps of the simulations
- A complete research paper (PDF) in the `paper/` directory
- A presentation (PDF) in the `presentation/` directory

## Paper and Presentation
### Building the Paper
1. To build paper, change directories into `paper/`
```
cd paper/
```
Then compile
```
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```
### Building the Presentation
2. To build presentation, change directories into `presentation/`
```
cd presentation/
```
1. Place all files in your working directory (or install to a local LaTeX path).
2. In your LaTeX document, load the theme with:

```latex
\documentclass{beamer}
\usetheme{CleanEasy}
```

3. (Optional) Include modular configurations:

```latex
\input{configs/configs.tex}
\input{configs/title_page.tex}
```

4. Compile with `pdflatex` or `latexmk`:

```bash
pdflatex CleanEasy.tex
```

---

## Contact
If you have any questions or suggestions, feel free to contact me at yuegelicay@gmail.com.

### Useful Resources
- [Henke X-ray database](https://henke.lbl.gov/) - Reference for X-ray transmission data


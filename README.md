# Kuang 2008 Linearized Model

[toc]

## Overview

This repository contains a Julia-based implementation of the linearized amtospheric model described in **Kuang (2008)**. The project is designed to simulate convective coupled waves with radiative feedback.

The framework couples a high-performance Julia simulation engine with a Python-based post-processing suite to calculate and visualize the results of the model.

## Repository Structure

The codebase is organized into simulation logic, configuration, experimental drivers, and analysis tools.

### Core Simulation (`src/`)
* `src/Kuang2008.jl`: The main entry point for the module.
* `src/LinearModel.jl`: Defining the operator matrix for the linearized model, including two-way and one-way modes.
* `src/params.jl`: Defining the physical and empirical parameters of the model, including four sets of parameters.
* `src/Diagnostics.jl`: Defining the diagnostic functions for the model.
* `src/Simulation.jl`: Defining the integration functions for the model.

### Configuration (`config/`)
* `configs/*.toml`:  TOML files defining experiments, and there are four regimes:
  * `no_rad.toml`: No radiative feedback.
  * `qt_rad.toml`: Radiative feedback on moisture and temperature.
  * `qt_cld_rad.toml`: Radiative feedback on moisture, temperature, and cloud.

### Experimental Drivers (`experiments/`)
* `experiments/run_case.jl`: Executing simulations.

### Data & Input (`data/`)
* `data/domain.h5`: Domain settings.
* `data/vertical_mode.h5`: Vertical modes.
* `data/inv_mat.h5`: Wavenumber samples and inverse Fourier matrix.
* `data/background.h5`: Background thermodynamic field.

### Post-Processing (`post/`)
* `post/kuang_post/make_all.py`: A batch processing script to generate all diagnostics for a run.
* `post/kuang_post/lib/`: Python modules for post-processing.
  * `post/kuang_post/lib/Diagnostics.py`: Calculates growth rates and phase speeds.
  * `post/kuang_post/lib/Plot.py`: Defining the plot functions for the model.
  * `post/kuang_post/lib/Reconstruct.py`: Defining the reconstruction functions for the model.

### Automation (`bin/ & Makefile`)
* `bin/run_case.sh`: Shell script to orchestrate the Julia simulation.

* `bin/postprocess.sh`: Shell script to trigger the Python analysis pipeline.

* `Makefile`: Facilitates building and workflow automation.

## Requirements

* Julia environment (version >= 1.8.5)
    Packages are managed via `Project.toml` and `Manifest.toml`
* Python environment (version >= 3.10.4)
    The virtual environment is saved in `post/`, and the Python dependencies are defined in `post/requirements.txt`
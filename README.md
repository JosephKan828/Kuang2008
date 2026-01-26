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

## How to Run

### Before execution

Before jump into the experiments, there are several things to check:

* Write a TOML file for you experiment in `configs/`.

  Within the directory, there is a `template.toml` file as an example. Please follow the instructure and comments to design your experiment.

* Ensure there is a parameter set corresponding to your experiment in `src/params.jl`:

  Since the parameter is interchangable for this model, there is a `struct` object for you to fill in all the necessary parameters. If you need additional parameters, you can append them to the struct.

* Check the domain settings in `data/domain.h5`:

  There is a pre-generated file for domain information, such as zonal and vertical coordinate. If you need to define new domain, background, vertical modes, or inverse matrix, please modify `Background.jl` for your own experiment.

### Execution

After you have check all the requirements, please run the following command:

* For running a single case (take `no_rad` as an example):

```bash
make run CASE=no_rad RAD_SCALE=0.001
```

In the above command, `no_rad` is the case name, and `0.001` is the radiative feedback scaling factor. It is worthwhile to notice that: **For `no_rad`** case, the radiative scaling is meaningless.

* For post processing a single case (take `no_rad` as an example):

```bash
make post CASE=no_rad
```

### Version Log

* v0.1.0: Initial release (2026/01/25)

### Contact Information

* Author: **Yu-Chuan Kan**

  Lab. of Chaos and Predictability,
  Department of Atmospheric Sciences,
  National Taiwan University,
  Taipei, Taiwan

* Email: [r14229003@ntu.edu.tw](mailto:r14229003@ntu.edu.tw) / [josephyck0828@gmail.com](mailto:josephyck0828@gmail.com) \

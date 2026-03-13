# When can negative control exposures be used for bias detection?

This repository contains the code for the simulations in the paper: "When can negative control exposures be used for bias detection?"

## Repository structure

This repository contains the R code used for the simulation studies and figures in  
**"When can negative control exposures be used for bias detection?"**

The intended entry point for using the code is the `scripts/` folder. These scripts run the main simulation scenarios by calling reusable functions stored in `R/`, saving outputs to `results/`, and producing figures from `figures/`.

### Main folders

- `scripts/`  
  Main entry point for the project. These scripts run the simulation scenarios analyzed in the paper.

- `R/`  
  Reusable functions for:
  - data generation
  - model fitting
  - parameter-grid creation
  - simulation execution

- `results/`  
  Saved simulation outputs, organized by scenario.

- `figures/`  
  Figure-making scripts and exported figure files.

### Workflow

A typical workflow is:

1. Run a scenario script from `scripts/`
2. The script sources functions from `R/`
3. Simulation outputs are saved to `results/`
4. Figure scripts in `figures/` use those results to create the manuscript figures

### Notes on naming

- `NoV` = scenario without the additional common cause `V`
- `V_D1` and `V_D2` = two scenarios with `V`
- `HighRho` = high-correlation variant
- `VaryA1Rho` = variant that changes how exposure–NCE correlation is generated
- `WithX` = extension with measured confounders `X`

### Acknowledgment

ChatGPT (OpenAI) was used to assist in drafting this README text.
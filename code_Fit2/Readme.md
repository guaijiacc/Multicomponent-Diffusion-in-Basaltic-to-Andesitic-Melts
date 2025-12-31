## Overview

This program uses the **BFGS optimization method** to fit multicomponent diffusion data, assuming a **single eigenvector matrix [Q]** and **seven sets of eigenvalues**.

To account for uncertainties in the determination of **interface positions** and **initial conditions**, the following parameters are treated as additional fitting variables:

- Interface position **x₀** (one x₀ for each diffusion couple experiment)
- Initial condition parameters: **w̄**, **Δw**, **s̄**, and **Δs**

---

## Diffusion Dataset

The diffusion data used for fitting include:

- **12 andesite–andesite diffusion couple experiments** and **1 basalt–andesite diffusion couple experiment** conducted in this study
- **28 basaltic diffusion couple experiments** from  
  Guo and Zhang (2018, 2020) and Bai and Zhang (2025)

The diffusion data are divided into **three categories based on experimental temperature**.  
For each temperature, two Excel files are provided:

- **Measured diffusion profiles**
- **1σ errors on measured oxide concentrations**

For example, at **1350 °C**:
- **`1350_RealData_8Comp.xlsx`** contains all measured diffusion profiles.
- **`1350_RealError_8Comp.xlsx`** contains the 1σ errors on measured oxide concentrations for all experiments.

---

## Regularization Parameters

The Excel file **`Penalty_factor.xlsx`** contains the **regularization hyperparameters L** assigned to each experiment.  
These parameters are used to penalize large absolute values of **s̄** and **Δs**.

- Unit of **L**: **5 × 10⁹**

---

## How to Use the Program

This folder contains the following MATLAB files:

1. **`BFGS_main_parallel.m`** — main program  
2. **`calculate_J_par.m`** — subroutine  
3. **`calculate_J7_par.m`** — subroutine  
4. **`con.m`** — subroutine  
5. **`prediction.m`** — subroutine  

The main program, **`BFGS_main_parallel.m`**, reads the initial guesses for:
- natural logarithm of eigenvalues,
- eigenvector matrix,
- interface positions,
- initial conditions,

from the Excel file **`Initial.xlsx`**.

---

## Parameter Definitions and Units

- **Eigenvalues**: μm²/s  
- **Eigenvector matrix elements**: dimensionless  
- **Interface positions (x₀)**: μm  
- **w̄**, **Δw**: wt%  
- **s̄**, **Δs**: wt%/μm  

---

## Output and Runtime

After approximately **30 minutes to 1 hour of runtime**, the optimized parameters and their **1σ uncertainties** can be accessed via the following variables:

- **Eigenvalues**: `trans_beta`, `error_trans_beta`
- **Interface positions**: `X0`, `error_X0`
- **Initial conditions**: `boundary`, `error_boundary`

The variable **`boundary`** is a **2D array**, where **every four columns** correspond to the optimized values of:

- w̄
- Δw
- s̄
- Δs  

for a single diffusion couple experiment.

---

## Output Files

The program generates an Excel file named:

- **`fitted_BA_AND.xlsx`**

which contains all optimized fitting parameters.

For reference, the folder also includes:

- **`fitted_BA_AND_from_author.xlsx`**

This file contains the optimized parameters generated on the authors’ computer.  
Users may compare their results with this file to verify consistency.

---

## Notes

For additional implementation details, numerical methods, and assumptions, users are encouraged to consult the **comments and notes within each MATLAB code file**.

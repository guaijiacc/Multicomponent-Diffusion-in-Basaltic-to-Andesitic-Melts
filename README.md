# Multicomponent Diffusion in Basaltic to Andesitic Melts  
**Toward Temperature- and Composition-Independent Diffusion Eigenvectors**

## Overview

This repository contains **open-access data files and MATLAB codes** associated with the study:

**Bai, B., Zhang, Y. (2026).**  
*Multicomponent diffusion in basaltic to andesitic melts: Toward a temperature- and composition-independent eigenvector matrix.*  
Submitted to **Geochimica et Cosmochimica Acta**.

The full archival version of this dataset is hosted on the  
**University of Michigan Deep Blue Data Repository**.  
**Deep Blue DOI / URL:** *(https://doi.org/10.7302/3dwy-7d91)*

---

## Scientific Background

Multicomponent diffusion in natural silicate melts controls many igneous processes, including magma mixing and crystal growth. In an **N-component** melt, diffusion is described by an **(N−1) × (N−1) diffusion matrix [D]**, whose:

- **Eigenvectors** define **N−1 eigen-components** that diffuse independently
- **Eigenvalues** represent the diffusion coefficients of these eigen-components

Previous studies demonstrated that diffusion eigenvectors in **8-component basaltic melts** (SiO₂–TiO₂–Al₂O₃–FeO–MgO–CaO–Na₂O–K₂O) are approximately **temperature independent** (Guo & Zhang, 2018, 2020; Bai & Zhang, 2025).

In this study, we extend the investigation to **andesitic compositions** to test whether diffusion eigenvectors remain invariant with **temperature**. We further explore the **composition** independence of diffusion eigenvectors in basaltic to andesitic melts.

---

## Experimental Dataset

A total of **13 diffusion couple experiments** were conducted at **1250 °C, 1350 °C, and 1500 °C**, including:

- **12 andesite–andesite diffusion couples**
- **1 basalt–andesite diffusion couple**

Diffusion profiles from the twelve andesite–andesite experiments can be **simultaneously fitted using a single eigenvector matrix**, supporting **temperature invariance** in andesitic melts.

The eigenvector matrix obtained for andesite is **similar to that in basalt**, suggesting possible **compositional invariance** across the basaltic–andesitic range.

To further test this hypothesis, diffusion profiles from:
- this study (andesitic and basalt–andesite couples), and
- **28 basaltic diffusion couples** ( 27 from Guo & Zhang, 2018, 2020 and 1 from Bai & Zhang, 2025)

are **simultaneously fitted using a single eigenvector matrix**.  
All diffusion profiles are well reproduced, supporting the hypothesis that diffusion eigenvectors are **both temperature- and composition-independent** from basaltic to andesitic melts.

---

## Methodology

1. **Diffusion eigenvectors and eigenvalues** are obtained using the **BFGS optimization method** implemented in MATLAB.

2. The fitting dataset includes:
   - 12 andesite–andesite diffusion couples (this study)
   - 1 basalt–andesite diffusion couple (this study)
   - 28 basaltic diffusion couples from previous studies

---

## Repository Contents

### Data Files

1. **`MultiComponentDif_calculator_v2.0.xlsx`**  
   An updated multicomponent diffusion calculator for computing diffusion profiles in **basaltic and andesitic melts**.

2. **`Bai_and_Zhang_2026_data_in_figures.xlsx`**  
   Data used to generate all figures presented in the associated manuscript.

---

### MATLAB Codes

3. **`code_Fit2.zip`**  
   MATLAB program (`BFGS_main_parallel.m`) for determining:
   - a single eigenvector matrix
   - seven sets of eigenvalues  

   The package includes all subroutines and diffusion data from **41 diffusion couple experiments** used in the fitting.  
   See the included `Readme` file for additional details.

---

## Requirements

- MATLAB (recent versions recommended)
- Optimization Toolbox (recommended)

---

## Related Publication

- **Bai, B., Zhang, Y. (2026).**  
  *Multicomponent diffusion in basaltic to andesitic melts: Toward a temperature- and composition-independent eigenvector matrix.*  
  (In progress)

---

## How to Cite

### Data citation (Deep Blue)

> Bai, B., Zhang, Y. (2026).  
> **Data of the paper "Multicomponent diffusion in basaltic to andesitic melts: Toward a temperature- and composition-independent eigenvector matrix"**  
> \[Data set\]. University of Michigan – Deep Blue Data.

**Deep Blue DOI / URL:** *(https://doi.org/10.7302/3dwy-7d91.)*

---

## Contact

**Bobo Bai**  
University of Michigan  
Email: bbai@umich.edu

# Py-RePS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%% Py-Reps - v. Feb 2025 %%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
Py-RePS is a python equivalent of RePS code(https://github.com/matteozennaro/reps)  
Latest version on https://github.com/vikhyat108/Py-RePS  


This repository contains two folders, `1-Fluid` and `2-Fluid`, which correspond to two different cosmological models:

1. **1-Fluid (Only CDM)**: A cosmology with total matter density consisting only of Cold Dark Matter (CDM) particles.
2. **2-Fluid (Massive Neutrinos + CDM)**: A cosmology that includes both CDM and massive neutrino particles.

## 1-Fluid (Only CDM)
The `1-Fluid` folder contains a Python script that generates the rescaled:
- Growth rates
- Transfer functions
- Power spectrum P(k)

These serve as the initial conditions for N-body simulations. Additionally, a parameter file allows users to modify cosmological parameters for their model.

**Note:** The `FF_GG` folder is not used in the 1-Fluid model; it is included for consistency purposes only.

## 2-Fluid (Massive Neutrinos + CDM Cosmology)
The `2-Fluid` folder also contains a Python script that performs the same tasks as in the `1-Fluid` model, generating:
- Growth rates
- Transfer functions
- Power spectrum P(k)

Users can modify the cosmology using the provided parameter file.

The `FF_GG` folder is used in this case to compute an integral required for the 2-Fluid model. The scripts for generating this file are inside the folder. However, users do not need to regenerate this file repeatedly—it can be generated once and used indefinitely, as it does not need to be updated for different cosmological models.

## Running the Code
To execute the script, use the following command:

```sh
python3 py-reps_1F.py
```

After execution, the output files will be saved in the `1F_ICs` folder.

---
For further modifications and parameter adjustments, refer to the parameter file within each respective folder.



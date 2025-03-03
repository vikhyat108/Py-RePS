# Py-RePS  
***********************************  
                ********* Py-Reps - v. Feb 2025 *********  
***********************************
  
Py-RePS is a python equivalent of RePS code(https://github.com/matteozennaro/reps)  
Latest version on https://github.com/vikhyat108/Py-RePS  

## Dependencies  
- Python package for CLASS i.e. "classy v3.2.1" (Can be installed using "pip install classy")  
- scipy 1.10.0  
- numpy 1.24.1  

*****************************************

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
If you have installed "classy v3.2.1" then you can directly use the following command:

```sh
python3 py-reps_1F.py  ####  For 1-Fluid case

python3 py-reps_2F.py  ####  For 2-Fluid case
```
  
  
Otherwise if you are having trouble with classy, you can choose to run CLASS(https://github.com/lesgourg/class_public version---2.10) code manually, for doing this do the following:-  
 - Copy CLASS code's folder inside the Py-RePS folder, where you have 1-Fluid and 2-Fluid folders (alternatively you can specify the path to your CLASS folder inside the py-reps_1F_class.py and py-reps_2F_class.py files respectively).
 - Rename the CLASS code's folder as "CLASS".
 - Now make two folders named "output_1F" and "output_2F" inside CLASS folder.  
 - Finally, copy the files named "parameters_1F_class.ini" and "parameters_2F_class.ini" from 1-Fluid and 2-Fluid folders respectively to the CLASS folder.

After doing this now run the command:

```sh
python3 py-reps_1F_class.py  ####  For 1-Fluid case

python3 py-reps_2F_class.py  ####  For 2-Fluid case

```
 



After execution, the output files will be saved in the `1F_ICs` and `2F_ICs` folders respectively.

---
For further modifications and parameter adjustments, refer to the parameter file within each respective folder.



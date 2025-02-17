# **Assignment 1: Elasto-Plasticity**  

## **Overview**  
This repository contains a Python implementation of the **Elasto-Plasticity models** based on **Kinematic** and **Isotropic** hardening.  

The project includes:  
- `elasto_plasticity.py` → The main containing modules related to kinematic and isotropic hardening. It also includes a plotting function to plot stress vs strain plots.
- `test_elasto_plasticity.py` → Unit tests using `pytest`  
- `tutorial.ipynb` → A Jupyter Notebook with five examples for demonstrating the functionality of the main module.

For information about Kinematic hardening and Isotropic hardening please visit this [website](https://www.fidelisfea.com/post/isotropic-kinematic-or-mixed-mode-which-hardening-model-for-your-abaqus-fea-analysis).

---

## **Installation & Requirements**  
Ensure you have **Python 3.12** installed along with the following dependencies:  

```bash
pip install numpy matplotlib pytest
```

---
**Please ensure that you run the following commands in the terminal after downloading the repository (Please ensure that repository is not in the downloads folder and their relative locations are not changed.)**\
To install the package, first create a virtual environment:
```bash
conda create --name elasto-plasticity-env python=3.12
```
Once the environment has been created, activate it:
```bash
conda activate elasto-plasticity-env
```
Double check that python is version 3.12 in the environment:
```bash
python --version
```
Please ensure that pip is using the most up to date version of setuptools:
```bash
pip install --upgrade pip setuptools wheel
```
Create the editable install of the Newton Solver code (note: you must be in the correct directory, i.e. where all of the files of the repository are.)
```bash
pip install -e .
```
Test the code is working with pytest
```bash
pytest -v --cov=elastoplasticity --cov-report term-missing
```

## **Usage**  
### **Importing the Module**  
To use the hardening models, first import the main module into your script or notebook:  

```python
from elastoplasticity import elasto_plasticity as ep
import numpy as np
from pathlib import Path
from typing import Callable
```

### **Initializing the different classes (Kinematic Hardening and Isotropic Hardening)**  
The classes can be initialized as:  

```python
'''
Kinematic Hardening
'''
ep.ElastoPlasticKinematic(E, Et, sigma_y)
```
```python
'''
Isotropic Hardening
'''
ep.ElastoPlasticIsotropic(E, Et, sigma_y)
```

Both classes need three input parameters for initialization.
#### **Parameters:**  
| Parameter   | Description |
|------------|-------------|
| `E`        | Young's modulus (Positive number) |
| `Et`        | Tangent modulus (Positive number smaller than Young's modulus) |
| `sigma_y`  | Yield stress (Positive number) |

**Note :** In case tangent modulus is not known but the hardening modulus is known then user needs to calculate it first.
#### **Example Initialization of both classes**  
```python
kin_model = ep.ElastoPlasticKinematic(E=1000, Et=100, sigma_y=10)
iso_model = ep.ElastoPlasticIsotropic(E=1000, Et=100, sigma_y=10)
```
### **Updating different parameters for the given strain history**

Once the classes have been initialized, the stress history for the given history can be generated by calling the update_steps function. For the above mentioned example, one can generate stresses for a given strain field as:

```python
strain_field = # Insert a NumPy array here

strain_increments = np.diff(strain_field, prepend=0)

# Step by step increase in strain
for delta_eps in strain_increments:
    kin_model.update_step(delta_eps)
    iso_model.update_step(delta_eps)
```
**Note :** Remember that update_step function takes increments in strain as an input not the strain field itself.

For understanding how the algorithm works for both classes, please visit [Piazza](https://cdn-uploads.piazza.com/paste/k7mp12pmibt1h5/b5b1c1452bb2d319f8da45dc127c5bf8422307ff85d4ae549bb49639fc9a64bd/Copy_of_Lecture_6_Blank.pdf).

### Plotting the stress-strain field
A plotting function has also been included in the main module. It can be called as follows for the above mentioned example:
```python
plot_stress_strain(kin_model, iso_model)
```
This function plots strain history to visualize the strain input. It also takes two parameters as an input which are associated with the stress and strain history of kinematic hardening model and isotropic hardening model.

---

## **Tutorials & Examples**  
To run the tutorial, user needs to install jupyter:
```bash
pip install jupyter
```
Then change the directory to tutorials folder:
```bash
cd tutorials/
```
Run the tutorial file:
```bash
jupyter notebook tutorial.ipynb
```
The **`tutorial.ipynb`** file provides five detailed examples to demonstrate the functionality of main module:  
#### Example 1: Lecture example
This example demonstrates the functionality of main module by plotting the stress-strain plot for the strain history that was solved during the [Lecture 6](https://cdn-uploads.piazza.com/paste/k7mp12pmibt1h5/b5b1c1452bb2d319f8da45dc127c5bf8422307ff85d4ae549bb49639fc9a64bd/Copy_of_Lecture_6_Blank.pdf). User can verify the results by comparing the generated plots to the lecture notes.
#### Example 2: Negative strains of that in example one.
This example demonstrates the functionality of the main module even if the system starts with compression.
#### Example 3: Sinusoidal strain
This example plots the stress-strain plots for the strain history input of the form $0.01\sin{x}$ for $x \in [0,4\pi]$.
The yield stress value has been altered from the previous examples so that yielding can take place easily.
#### Example 4: Modulated Sinusoidal strain
This example plots the stress-strain plots for the strain history input of the form $0.01x\sin{x}$ for $x \in [0, 8\pi]$. Yield stress has also been altered for this example so that the material doesn't yield at the beginning.
#### Example 5: Exponential loading and unloading
The strain history is taken of the form of piecewise exponential functions for loading and unloading.

---

## **Error Handling & Warnings**  
The module raises errors in the following cases:  
- If Yield stress is negative (`ValueError: Yield stress should be positive.`)  
- If Young's modulus is not larger than Tangent modulus (`ValueError: Invalid input: Et {Et} is either greater than or equal to E {E}.`)

## Final remarks for the user
User is encouraged to play with different parameters in *test_elasto_plasticity.py* and *tutorial.ipynb* notebook to check how the script *elasto_plasticity.py* works under different scenarios.

'''
This file contains kinematic hardening and isotropic hardening models implemented as separate classes to study the elasto-plastic
behavior of a material.
'''
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Callable, Union

'''
Hello world function to check the proper import of the module.
'''
def hello_world():
    return "hello world!"
    
'''
Check for Young's modulus being larger than the tangent modulus.
'''
def check_E_larger_than_Et(E,Et):
    if E > Et:
        return True
    else:
        raise ValueError(f"Invalid input: Et {Et} is either greater than or equal to E {E}.")
        
'''
Check if the yield stress is positive or not
'''
def check_valid_sigma_y(sigma_y):
    if sigma_y >= 0:
        return True
    else:
        raise ValueError(f"Yield stress should be positive.")

'''
Implemenation of Elasto-Plastic model with kinematic hardening
'''
class ElastoPlasticKinematic:
    '''
    Elasto-Plastic Model with Kinematic Hardening
    E: Young's Modulus
    H: Hardening Modulus
    sigma_y: Yield Strength
    eta: Hardening rate (plastic multiplier)
    '''
    def __init__(self, E, Et, sigma_y):
        self.E = E  # Young's Modulus
        self.Et = Et  # Tangent modulus
        check_E_larger_than_Et(self.E, self.Et)
        self.H = self.E * self.Et / (self.E - self.Et)  # Hardening modulus (based on Young's and Tangent moduli)
        self.sigma_y = sigma_y  # Initial yield strength
        check_valid_sigma_y(sigma_y)
        self.sigma_n = 0  # Current stress
        self.epsilon_n = 0  # Total strain
        self.epsilon_p_n = 0  # Plastic strain
        self.alpha_n = 0  # Kinematic hardening variable (back stress)
        self.eta = 0  # Hardening rate (plastic multiplier)
        self.strain_history = []  # Store strain values
        self.stress_history = []  # Store stress values
    '''
    Calculation of stress, yield stress, etc. with the changing strain field
    '''
    def update_step(self, delta_eps):
        # Trial stress (assumes stress would increase if no yielding occurs)
        sigma_trial = self.sigma_n + self.E * delta_eps  # Trial stress based on current strain increment

        # Define alpha_trial (back stress at trial state)
        alpha_trial = self.alpha_n

        # Define eta_trial as the difference between sigma_trial and alpha_trial
        eta_trial = sigma_trial - alpha_trial
        self.eta = eta_trial

        # Yield function: f_trial = |eta_trial| - sigma_y
        f_trial = abs(self.eta) - self.sigma_y  # Trial yield function based on eta_trial

        # Elastic regime
        if f_trial <= 0:
            self.sigma_n = sigma_trial
        # Plastic regime
        else:  
            delta_gamma = f_trial / (self.E + self.H)
            # Loading (positive increment in strain)
            if np.sign(sigma_trial) == np.sign(self.sigma_n):
                self.sigma_n = sigma_trial - np.sign(self.eta) * self.E * delta_gamma
                self.alpha_n = self.alpha_n + np.sign(self.eta) * self.H * delta_gamma
                self.epsilon_p_n = self.epsilon_p_n + delta_gamma
            # Unloading (negative increment in strain)
            else:
                self.sigma_n = sigma_trial - np.sign(self.eta) * self.E * delta_gamma
                self.alpha_n = self.alpha_n - np.sign(self.eta) * self.H * delta_gamma
                self.epsilon_p_n = self.epsilon_p_n - delta_gamma
        # Update the total strain
        self.epsilon_n += delta_eps
        # Store the strain and stress for future analysis (plotting or diagnostics)
        self.strain_history.append(self.epsilon_n)
        self.stress_history.append(self.sigma_n)
        # Return the current stress value
        return self.sigma_n
'''
Implemenation of Elasto-Plastic model with kinematic hardening
'''
class ElastoPlasticIsotropic:
    '''
    Elasto-Plastic Model with Isotropic Hardening
    E: Young's Modulus
    H: Hardening Modulus
    sigma_y: Yield Strength
    '''
    def __init__(self, E, Et, sigma_y):
        self.E = E      # Young's modulus
        self.Et = Et    # Tangent modulus
        check_E_larger_than_Et(self.E, self.Et)
        self.H = self.E * self.Et / (self.E - self.Et)  # Hardening modulus (based on Young's and Tangent moduli)
        self.sigma_y = sigma_y  # Initial yield strength
        check_valid_sigma_y(sigma_y)
        self.epsilon_n = 0  # Total strain
        self.epsilon_p_n = 0  # Plastic strain
        self.sigma_n = 0  # Current stress
        self.strain_history = []  # Store strain values
        self.stress_history = []  # Store stress values
    '''
    Calculation of stress, yield stress, etc. with the changing strain field
    '''
    def update_step(self, delta_eps):
        # Trial stress (assumes stress would increase if no yielding occurs)
        sigma_trial = self.sigma_n + self.E * delta_eps  # Trial stress based on current strain increment

        # Yield function: f_trial = |sigma_trial| - sigma_y
        f_trial = abs(sigma_trial) - self.sigma_y  # Absolute value for tension and compression
        
        # Elastic regime
        if f_trial <= 0:
            self.sigma_n = sigma_trial
        # Plastic regime    
        else:
            # Loading (Positive increment in strain)
            if np.sign(sigma_trial) == np.sign(self.sigma_n):  
                delta_gamma = f_trial / (self.E + self.H)  # Plastic multiplier
                self.sigma_n = sigma_trial - np.sign(sigma_trial) * self.E * delta_gamma  # Corrected stress
                self.epsilon_p_n += delta_gamma  # Increase plastic strain
                # Update yield stress due to hardening
                self.sigma_y += self.H * delta_gamma  # Isotropic hardening rule for loading
            # Unloading (Negative increment in strain)
            else:
                delta_gamma = f_trial / (self.E + self.H)  # Plastic multiplier during unloading
                self.sigma_n = sigma_trial + np.sign(sigma_trial) * self.E * delta_gamma  # Corrected stress
                self.epsilon_p_n -= delta_gamma  # Decrease plastic strain during unloading
                # Update yield stress due to hardening (but in the opposite direction for unloading)
                self.sigma_y -= self.H * delta_gamma  # Isotropic hardening rule for unloading

        # Update the total strain
        self.epsilon_n += delta_eps

        # Store the strain and stress for future analysis (plotting or diagnostics)
        self.strain_history.append(self.epsilon_n)
        self.stress_history.append(self.sigma_n)
        # Return the current stress value
        return self.sigma_n

def plot_stress_strain(kinematic_model, isotropic_model):
    '''
    Plot strain vs load stress and stress vs strain curves for both kinematic and isotropic hardening models.

    Parameters:
    kinematic_model (ElastoPlasticKinematic): Instance of kinematic hardening model
    isotropic_model (ElastoPlasticIsotropic): Instance of isotropic hardening model
    '''
    # Create the figure and axes
    fig, axs = plt.subplots(1, 3, figsize=(15, 7))

    # Strain plot
    axs[0].plot(isotropic_model.strain_history, 
                label='Strain history', color='red', linewidth=2.5)
    axs[0].set_title('Strain History', fontsize=25, fontweight='bold')
    axs[0].set_xlabel(r'Load Step', fontsize=25, fontweight='bold')
    axs[0].set_ylabel(r'Strain', fontsize=25, fontweight='bold')
    axs[0].grid(True)
    axs[0].tick_params(axis='both', which='major', labelsize=15)
    
    # Kinematic Hardening Plot
    axs[1].plot(kinematic_model.strain_history, kinematic_model.stress_history, 
                label='Kinematic Hardening', color='green', linewidth=2.5)
    axs[1].set_title('Kinematic Hardening', fontsize=25, fontweight='bold')
    axs[1].set_xlabel(r'Strain', fontsize=25, fontweight='bold')
    axs[1].set_ylabel(r'Stress', fontsize=25, fontweight='bold')
    axs[1].grid(True)
    axs[1].tick_params(axis='both', which='major', labelsize=15)
    
    # Isotropic Hardening Plot
    axs[2].plot(isotropic_model.strain_history, isotropic_model.stress_history, 
                label='Isotropic Hardening', color='blue', linewidth=2.5)
    axs[2].set_title('Isotropic Hardening', fontsize=25, fontweight='bold')
    axs[2].set_xlabel(r'Strain', fontsize=25, fontweight='bold')
    axs[2].grid(True)
    axs[2].tick_params(axis='both', which='major', labelsize=15)
    
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

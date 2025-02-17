'''
Test file for elastoplasticity modules
'''
from elastoplasticity import elasto_plasticity as ep
from io import BytesIO
import numpy as np
from pathlib import Path
import pytest
import re
import matplotlib.pyplot as plt

'''
Test 1: Checks the import of elasoplasticity module.
'''
def test_hello_world():
    known = "hello world!"
    assert known == ep.hello_world()

'''
Test 2: Tests if the function that compares the value of E and Et is functioning properly.
'''
def test_check_E_larger_than_Et():
    E = 100
    Et = 10
    assert ep.check_E_larger_than_Et(E, Et) == True

'''
Test 3: Tests if the code returns a ValueError when unreasonable values of E and Et are provided.
'''
def test_check_E_not_larger_than_Et():
    E = 10
    Et = 100
    with pytest.raises(ValueError, match=f"Invalid input: Et {Et} is either greater than or equal to E {E}."):
        ep.check_E_larger_than_Et(E, Et)

'''
Test 4: Tests if function is properly getting the reasonable values for yield stress.
'''
def test_valid_sigma_y():
    sigma_y = 1
    assert ep.check_valid_sigma_y(sigma_y) == True

'''
Test 5: Tests if the code returns a ValueError when unreasonable value of yield stress is provided.
'''
def test_invalid_sigma_y():
    sigma_y = -1
    with pytest.raises(ValueError, match=f"Yield stress should be positive."):
        ep.check_valid_sigma_y(sigma_y)

'''
Creating a dummy model to check the other functionalities of the main module.
'''
@pytest.fixture
def setup_models():
    # Create sample parameters for testing
    E = 210e9  # Young's modulus
    Et = 1e9   # Tangent modulus
    sigma_y = 250e6  # Yield stress
    
    # Initialize models
    kinematic_model = ep.ElastoPlasticKinematic(E, Et, sigma_y)
    isotropic_model = ep.ElastoPlasticIsotropic(E, Et, sigma_y)
    
    return kinematic_model, isotropic_model

'''
Test 6: Tests if the initial values of various parameters set to zero at the time of initialization.
'''
def test_initial_conditions(setup_models):
    kinematic_model, isotropic_model = setup_models
    
    # Test that the initial conditions are set correctly
    assert kinematic_model.epsilon_p_n == 0
    assert kinematic_model.alpha_n == 0
    assert isotropic_model.epsilon_p_n == 0
    assert isotropic_model.epsilon_n == 0
'''
Test 7: Tests if the values are updating properly or not for kinematic hardening during loading.
'''
def test_update_step_kinematic_loading(setup_models):
    kinematic_model, _ = setup_models
    
    # Apply a positive strain increment and test for loading behavior
    kinematic_model.update_step(0.001)
    
    # Check that stress has been updated
    assert len(kinematic_model.stress_history) > 0
    assert kinematic_model.stress_history[-1] > 0  # Stress should be positive after loading

'''
Test 8: Tests if the values are updating properly or not for isotropic hardening during loading.
'''
def test_update_step_isotropic_loading(setup_models):
    _, isotropic_model = setup_models
    
    # Apply a positive strain increment and test for isotropic loading behavior
    isotropic_model.update_step(0.001)
    
    # Check that stress has been updated
    assert len(isotropic_model.stress_history) > 0
    assert isotropic_model.stress_history[-1] > 0  # Stress should be positive after loading

'''
Test 9: Tests if the values are updating properly or not for kinematic hardening during unloading.
'''
def test_update_step_kinematic_unloading(setup_models):
    kinematic_model, _ = setup_models
    
    # Apply a positive strain increment (loading)
    kinematic_model.update_step(0.001)
    
    # Now apply a negative strain increment (unloading)
    kinematic_model.update_step(-0.001)
    
    # Check that stress has decreased
    assert kinematic_model.stress_history[-1] < kinematic_model.stress_history[-2]  # Stress should decrease after unloading

'''
Test 10: Tests if the values are updating properly or not for isotropic hardening during unloading.
'''
def test_update_step_isotropic_unloading(setup_models):
    _, isotropic_model = setup_models
    
    # Apply a positive strain increment (loading)
    isotropic_model.update_step(0.001)
    
    # Now apply a negative strain increment (unloading)
    isotropic_model.update_step(-0.001)
    
    # Check that stress has decreased
    assert isotropic_model.stress_history[-1] < isotropic_model.stress_history[-2]  # Stress should decrease after unloading

'''
Test 11: Tests if the plotting function is behaving properly or not.
'''
def test_plot_stress_strain(setup_models):
    kinematic_model, isotropic_model = setup_models

    # Apply some strain steps to both models
    for eps in [0.001, 0.002, 0.003]:
        kinematic_model.update_step(eps)
        isotropic_model.update_step(eps)

    # Test the plotting function by capturing the plot in a BytesIO object
    plt.ioff()  # Turn off interactive plotting to test without displaying

    # Create a plot in a BytesIO buffer
    buf = BytesIO()

    # Ensure the plot is created
    ep.plot_stress_strain(kinematic_model, isotropic_model)
    
    # Show the plot for debugging (temporary)
    plt.show()  # Temporarily check if the plot appears
    
    # Save the figure to the buffer
    plt.savefig(buf)
    buf.seek(0)  # Move to the beginning of the buffer

    # Check that the plot was created by verifying the buffer is not empty
    content = buf.read()
    assert content  # If the plot was generated, this should not raise an error

'''
Test 12: Tests if the plastic multipliers are functionung properly during unloading.
'''
def test_unloading_plastic_multiplier(setup_models):
    kinematic_model, isotropic_model = setup_models

    # Simulate loading phase (increasing strain)
    kinematic_model.update_step(0.01)  # Apply a small strain increment
    isotropic_model.update_step(0.01)

    # Check if the model has entered the plastic regime (should not be in elastic)
    assert kinematic_model.stress_history[-1] > kinematic_model.sigma_y, "The model didn't yield"
    
    # Simulate unloading phase (reversing the strain direction)
    kinematic_model.update_step(-0.01)  # Apply a negative strain increment
    isotropic_model.update_step(-0.01)

    # Ensure the unloading phase was correctly handled
    assert kinematic_model.stress_history[-1] < kinematic_model.sigma_y, "Unloading didn't occur correctly"

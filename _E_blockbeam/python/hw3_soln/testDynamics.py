import sys
import os
# Get the directory one level above
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add it to sys.path
sys.path.insert(0, parent_dir)
# use this file to test your dynamics file
import numpy as np
from blockBeamDynamics import blockBeamDynamics as dynamics
import blockbeamParam as P
import testCases as TC



# define test parameters
P.Ts = .01
plant = dynamics()
numTests = 0
numTestsPassed = 0
    
# test all test cases
for initial_state, input, true_state in TC.dynamics_test_cases:
    numTests += 1
    initial_state = np.array(initial_state)
    true_state = np.array(true_state)
    plant.state = initial_state
    plant.update(np.array(input))
    actual_state = plant.state
    error = np.linalg.norm(actual_state - true_state)
    if error < TC.errorTolerance:
        numTestsPassed += 1
        print(f'Test #{numTests} passed!')
        continue
    print(f'Test #{numTests} failed with an error magnitude of {error}')
    print(f"\tThe actual output was {actual_state.flatten()}")
    print(f"\tIt should've been: {true_state.flatten()}")

# display final results of the tests 
if numTestsPassed == numTests:
    print("\nExcellent work!! Your dynamics file has passed all of the tests!\n")
else:
    print(f"\nDarn, only {numTestsPassed}/{numTests} of the tests passed... Keep it up! You'll get it!\n")

import sys
import os

# Get the directory one level above
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add it to sys.path
sys.path.insert(0, parent_dir)

import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics

import matplotlib
matplotlib.use('tkagg')  # requires TkInter

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.02)
force = signalGenerator(amplitude=.5, frequency=1.0, y_offset=11.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        r = reference.square(t)
        u = force.sin(t)
        y = blockBeam.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockBeam.state)
    dataPlot.update(t, blockBeam.state, u, r) # Switched order for current dataPlotter -- Christian Hales, 9/23/2024
    plt.pause(0.0001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
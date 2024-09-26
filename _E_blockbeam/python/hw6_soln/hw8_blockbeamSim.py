import matplotlib.pyplot as plt
import blockBeamParam as P
from blockBeamDynamics import blockBeamDynamics
from ctrlPD import ctrlPD
from signalGenerator import signalGenerator
from blockBeamAnimation import blockBeamAnimation
from dataPlotter import dataPlotter

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics()
controller = ctrlPD()

#HW asked for reference input frequency of 0.01, but this is more interesting
reference = signalGenerator(amplitude=0.15, frequency=0.05, y_offset=0.25)
disturbance = signalGenerator(amplitude=0.0, frequency=0.0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockBeamAnimation()
t = P.t_start  # time starts at t_start
y = blockBeam.h()  # output of system at start of simulation
while t < P.t_end:  # main simulation loop

    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        r = reference.square(t)  # reference input
        d = disturbance.step(t)  # input disturbance
        n = 0.0  #noise.random(t)  # simulate sensor noise
        x = blockBeam.state
        u = controller.update(r, x)  # update controller
        y = blockBeam.update(u + d)  # propagate system
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockBeam.state)
    dataPlot.update(t, r, blockBeam.state, u)
    plt.pause(0.001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
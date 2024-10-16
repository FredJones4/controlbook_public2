# %% [markdown]
# # Midterm 1 Example Problems Solution Notebook
# 
# This notebook provides solutions to the pre-midterm questions using Python libraries (`numpy` and `sympy`), as allowed. 
# Each problem is solved in its respective section with appropriate explanations.
# 
# ## Question 1: Find the Kinetic Energy, Potential Energy, and EOM for the system in Figure 1.

# %%
import numpy as np
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy import sin, cos

# Define the symbols and variables
t = sp.symbols('t')  # time
m1, m2, mb, g, l = sp.symbols('m1 m2 mb g l')  # masses, gravity, length
z = dynamicsymbols('z')  # vertical displacement
theta = dynamicsymbols('theta')  # angular displacement
z_dot = z.diff(t)  # time derivative of z
theta_dot = theta.diff(t)  # time derivative of theta

# Kinetic Energy: Translational and Rotational
T_trans = (1/2) * (m1 + m2 + mb) * z_dot**2  # Translational KE
T_rot = (1/2) * (m1 * l**2 + m2 * l**2) * theta_dot**2  # Rotational KE (assuming simple rotation)
T = T_trans + T_rot

# Potential Energy: Gravitational
V = (m1 + m2 + mb) * g * z * sp.cos(theta)  # Gravitational potential

# Lagrangian
L = T - V

# Display Lagrangian
sp.pprint(L)

# %% [markdown]
# ### Equation of Motion (EOM)
# 
# Using the Euler-Lagrange equation to derive the EOMs for \( z \) and \( \theta \).

# %%
# Using the derive_Lagrangian function defined earlier
def derive_Lagrangian(L, q, q_dot):
    term_1 = (sp.diff(L, q_dot)).diff(t)
    term_2 = sp.diff(L, q)
    return term_1 - term_2

# Derive EOM for z and theta
EOM_z = derive_Lagrangian(L, z, z_dot)
EOM_theta = derive_Lagrangian(L, theta, theta_dot)

# Display EOMs
sp.pprint(EOM_z)
sp.pprint(EOM_theta)

# %% [markdown]
# ## Question 2: Generalized Equilibrium Forces for the System
# 
# Given the provided equations of motion:
# 
# \[
# (m_b + m_1 + m_2) \ddot{z} + \frac{\ell}{2}(m_1 - m_2) \cos \theta \ddot{\theta} = -b \dot{z} + \tau + \frac{\ell}{2} (m_1 - m_2) \dot{\theta}^2 \sin \theta
# \]
# \[
# \frac{\ell^2}{4} \left( \frac{m_b}{3} + m_1 + m_2 \right) \ddot{\theta} + \frac{\ell}{2} (m_1 - m_2) \ddot{z} \cos \theta = F \cos \theta + g(m_1 + m_2)
# \]

# %% [markdown]
# We will substitute the equilibrium conditions \( \ddot{z} = 0 \), \( \ddot{\theta} = 0 \), \( \dot{z} = 0 \), and \( \dot{\theta} = 0 \) into the given equations to find the generalized forces.

# %%
# Define the additional parameters
b, tau, F = sp.symbols('b tau F')

# Generalized Equilibrium Equations
eq1 = (mb + m1 + m2) * z.diff(t, t) + (l/2) * (m1 - m2) * cos(theta) * theta.diff(t, t) + b * z.diff(t) - tau
eq2 = (l**2 / 4) * (mb/3 + m1 + m2) * theta.diff(t, t) + (l/2) * (m1 - m2) * z.diff(t, t) * cos(theta) - F * cos(theta)

# Simplify at equilibrium (z_dot = 0, theta_dot = 0, z_ddot = 0, theta_ddot = 0)
equilibrium_eq1 = eq1.subs({z.diff(t): 0, z.diff(t, t): 0, theta.diff(t): 0, theta.diff(t, t): 0})
equilibrium_eq2 = eq2.subs({z.diff(t): 0, z.diff(t, t): 0, theta.diff(t): 0, theta.diff(t, t): 0})

# Display equilibrium forces
sp.pprint(equilibrium_eq1)
sp.pprint(equilibrium_eq2)

# %% [markdown]
# ## Question 3: Transfer Function
# 
# Given the differential equation:
# 
# \[
# \dddot{y} + 2\ddot{y} + 3\dot{y} + 4y = 5\dot{u} + 6u
# \]
# 
# We'll compute the transfer function by taking the Laplace transform.

# %%
# Define Laplace variable
s = sp.symbols('s')

# Define Y(s) and U(s) in the Laplace domain
Y = sp.Function('Y')(s)
U = sp.Function('U')(s)

# Transfer Function for the given DE
transfer_func = (5*s + 6) / (s**3 + 2*s**2 + 3*s + 4)

# Display transfer function
sp.pprint(transfer_func)

# %% [markdown]
# ## Question 4: Stability Analysis
# 
# Given the closed-loop transfer function:
# 
# \[
# G(s) = \frac{s - 7}{2s^2 - 4s + 54}
# \]
# 
# We'll analyze the denominator to classify stability.

# %%
# Define the closed-loop transfer function
G = (s - 7) / (2*s**2 - 4*s + 54)

# Find the characteristic equation (denominator)
char_eq = sp.denom(G)

# Find the roots of the characteristic equation
roots = sp.solve(char_eq, s)

# Display roots
sp.pprint(roots)

# %% [markdown]
# ## Question 5: Closed-Loop Characteristic Equation for PD Control
# 
# Given the system:
# 
# \[
# P(s) = \frac{8}{2s^2 + 4s - 6}
# \]
# 
# The closed-loop characteristic equation for PD control is found by solving \( 1 + P(s)C(s) = 0 \).

# %%
# Define P(s)
P = 8 / (2*s**2 + 4*s - 6)

# Assume a PD controller of the form C(s) = Kp + Kd*s
Kp, Kd = sp.symbols('Kp Kd')
C = Kp + Kd * s

# Closed-loop characteristic equation
closed_loop_eq = 1 + P * C

# Display characteristic equation
sp.pprint(sp.expand(closed_loop_eq))

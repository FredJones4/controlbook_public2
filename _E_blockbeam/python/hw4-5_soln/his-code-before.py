# %% [markdown]
# # E.3 Calculations

# %%
import sympy
from sympy import *
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

init_printing()

def dotprint(expr):
    display(Math(vlatex(expr)))

# %% [markdown]
# ### Define Symbols and Configuration Variables

# %%
# Time symbol
t = symbols('t')

# %%
# Generalized coordinates

# TODO: Change this!
z = symbols('z', cls=Function)
z = z(t)
theta = symbols('\\theta', cls=Function)
theta = theta(t)

# %%
# Derivatives of Generalized Coordinates

# TODO: Change this!
z_dot = z.diff(t)
theta_dot = theta.diff(t)

# %%
# Other symbols (mass, forces, damping coefficients, torques...)

# TODO: Change this!
m1, m2, F, el, g  = symbols('m_1, m_2, F, \\ell, g')

# %% [markdown]
# ### Kinetic Energy
# 
# #### Translational Kinetic Energy

# %%
# TODO: Change this!
K_T = 1/6*el**2*m2*theta_dot**2 + 1/2*m1*z**2*theta_dot**2 + 1/2*m1*z_dot**2

dotprint(K_T)

# %% [markdown]
# #### Rotational Kinetic Energy (ignore)

# %%
# TODO: Change this!
K_R = 0

dotprint(K_R)

# %% [markdown]
# #### Total Kinetic Energy

# %%
K = K_T #+ K_R

dotprint(K)

# %% [markdown]
# ### Potential Energy

# %%
# TODO: Change this!
P = m1*g*z*sin(theta) + m2*g*el*sin(theta)

dotprint(P)

# %% [markdown]
# ### Define the Lagrangian
# 
# For this class, the Lagrangian will in general be a scalar value equal to the difference between kinetic and potential energy.

# %%
L = K - P

dotprint(L)

# %% [markdown]
# ### Define the Forces
# 
# #### Generalized Forces
# 
# The forces will form a vector with the same length as the number of configuration variables. 
# 
# The first entry corresponds to the first configuration variable, the second entry to the second, etc.
# 
# Each entry should contain the generalized forces acting on the corresponding configuration variable.

# %%
# TODO: Change this!
tau = Matrix([0, F*el*cos(theta)])

dotprint(tau)

# %% [markdown]
# #### Non-conservative Forces (e.g. Damping)
# 
# Similar to the generalized forces, this will be a vector with the same length as the number of configuration variables.
# 
# Each entry should contain the nonconservative forces acting on the corresponding configuration variable.

# %%
# TODO: Change this!
B = Matrix([0, 0])
q_dot = Matrix([0])

Bqdot = B @ q_dot

dotprint(Bqdot)

# %%
F_total = tau + Bqdot

dotprint(F_total)

# %% [markdown]
# ### Compute the Lagrangian Derivatives / Partial Derivatives
# 
# The following function computes these two terms:
# 
# * $\frac{d}{dt} \left(\frac{\partial L(q, \dot{q})}{\partial \dot{q}} \right)$
# * $\frac{\partial L(q,\dot{q})}{\partial q}$
# 
# Things to note:
# 
# * `q` and `q_dot` should be Python tuples containing your configuration variables (e.g. `(x, theta)`)
# * The function returns a Python list. You probably want to wrap the output in a Sympy `Matrix()`. For example,
# 
# ```python
# out = derive_Lagrangian(...)
# out = Matrix(out)
# ```

# %%
def derive_Lagrangian(L, z, theta, z_dot, theta_dot):
    term_1 = (sympy.tensor.derive_by_array(L, z_dot)).diff(t)
    term_2 = sympy.tensor.derive_by_array(L, z)
    term_3 = (sympy.tensor.derive_by_array(L, theta_dot)).diff(t)
    term_4 = sympy.tensor.derive_by_array(L, theta)
    ans = Matrix([term_1 - term_2, term_3 - term_4])
    return ans

# %% [markdown]
# ### Get the final Euler Lagrange Equations of Motion
# 
# Remember that the Euler Lagrange Equations used in class are:
# 
# $$
# \underbrace{\frac{d}{dt} \left(\frac{\partial L(q, \dot{q})}{\partial \dot{q}} \right) - \frac{\partial L(q,\dot{q})}{\partial q}}_\text{Left Hand Side (LHS)} = \underbrace{\tau - B \dot{q}}_\text{Right Hand Side (RHS)}.
# $$
# 
# In Sympy, we use the `Eq(LHS, RHS)` class to get a symbolic equation $\text{LHS} = \text{RHS}$.
# 
# Both the left and right expressions should be Sympy vectors (i.e. Matrices with 1 column). 
# 
# The length of the left and right vectors must both be equal to the number of configuration variables.

# %%
# TODO: Change this!
LHS = derive_Lagrangian(L, z, theta, z_dot, theta_dot)
RHS = F_total

Euler_Lagrange = Eq(LHS, RHS)

dotprint(simplify(Euler_Lagrange))



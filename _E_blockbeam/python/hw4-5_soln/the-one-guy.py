# %% [markdown]
# # E.3 Calculations and Linearization

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
z = symbols('z', cls=Function)
z = z(t)
theta = symbols('\\theta', cls=Function)
theta = theta(t)

# %%
# Derivatives of Generalized Coordinates
z_dot = z.diff(t)
theta_dot = theta.diff(t)

# TODO: Add second derivatives for linearization
z_ddot = z.diff(t, 2)
theta_ddot = theta.diff(t, 2)

# %%
# Other symbols (mass, forces, damping coefficients, torques...)
m1, m2, F, el, g  = symbols('m_1, m_2, F, \\ell, g')

# %% [markdown]
# ### Kinetic Energy
# 
# #### Translational Kinetic Energy

# %%
K_T = 1/6*el**2*m2*theta_dot**2 + 1/2*m1*z**2*theta_dot**2 + 1/2*m1*z_dot**2

dotprint(K_T)

# %% [markdown]
# #### Rotational Kinetic Energy (ignore)

# %%
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
P = m1*g*z*sin(theta) + m2*g*el*sin(theta)

dotprint(P)

# %% [markdown]
# ### Define the Lagrangian

# %%
L = K - P

dotprint(L)

# %% [markdown]
# ### Define the Forces
# 
# #### Generalized Forces

# %%
tau = Matrix([0, F*el*cos(theta)])

dotprint(tau)

# %% [markdown]
# #### Non-conservative Forces (e.g. Damping)

# %%
B = Matrix([0, 0])
q_dot = Matrix([z_dot, theta_dot])  # TODO: Changed this to include both generalized coordinates

Bqdot = B @ q_dot

dotprint(Bqdot)

# %%
F_total = tau + Bqdot

dotprint(F_total)

# %% [markdown]
# ### Compute the Lagrangian Derivatives / Partial Derivatives

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

# %%
LHS = derive_Lagrangian(L, z, theta, z_dot, theta_dot)
RHS = F_total

Euler_Lagrange = Eq(LHS, RHS)

dotprint(simplify(Euler_Lagrange))

# TODO: Add linearization steps

# %% [markdown]
# ## Deriving Nonlinear State Space Equations

# %%
# Define the state vector
state = Matrix([z, theta, z_dot, theta_dot])
dotprint(state)

# %%
# Solve the Euler-Lagrange equations for z_ddot and theta_ddot
solved = solve(Euler_Lagrange, (z_ddot, theta_ddot))

# Define the state derivative vector
state_deriv = Matrix([
    z_dot,
    theta_dot,
    solved[z_ddot],
    solved[theta_ddot]
])
dotprint(state_deriv)

# %% [markdown]
# ## Linearization
# 
# First, we need to find an equilibrium point. We can do this by setting all derivatives to zero and solving.

# %%
equilibrium_equation = state_deriv.subs({z_dot: 0, theta_dot: 0, theta: 0})

eq_solve_dict = solve(equilibrium_equation, (z, F), simplify=True, dict=True)[0]
dotprint(eq_solve_dict)

# %% [markdown]
# We can see that at equilibrium, $z$ can be any value (which we'll call $z_e$), and $F$ depends on this $z_e$.

# %%
z_e = symbols('z_e')
u_eq = m1*g/el*z_e + m2*g/2
theta_eq = 0
z_dot_eq = 0
theta_dot_eq = 0

# %% [markdown]
# #### Define A, B Jacobians
# 
# We can use Sympy's `jacobian` function to find the jacobians of `f(x,u)`.
# 
# First we find $A = \frac{\partial f}{\partial x}$:

# %%
A = state_deriv.jacobian(state)
dotprint(A)

# %%
A_subs = {
    z: z_e,
    theta: theta_eq,
    z_dot: z_dot_eq,
    theta_dot: theta_dot_eq,
    F: u_eq
}

A_eq = A.subs(A_subs)
dotprint(A_eq)

# %% [markdown]
# Now we do a similar process to find $B = \frac{\partial f}{\partial u}$

# %%
B = state_deriv.jacobian(Matrix([F]))
dotprint(B)

# %%
B_subs = {
    z: z_e,
    theta: theta_eq,
    z_dot: z_dot_eq,
    theta_dot: theta_dot_eq,
    F: u_eq
}

B_eq = B.subs(B_subs)
dotprint(B_eq)

# %% [markdown]
# ### Transfer Function
# 
# We can also transform this to a transfer function if we define C and D matrices

# %%
C = Matrix([[0, 1, 0, 0], [0, 0, 0, 1.0]])
D = Matrix([[0], [0]])

s = symbols('s')
transfer_func = simplify(C * (s*eye(4) - A_eq).inv() * B_eq + D)
dotprint(transfer_func)

# %% [markdown]
# ### Simplifying Assumption
# 
# Now setting the $m_1g$ term equal to zero as described in the problem:

# %%
A_simplified = A_eq.subs(m1*g, 0)
B_simplified = B_eq.subs(m1*g, 0)

transfer_func_simplified = simplify(C * (s*eye(4) - A_simplified).inv() * B_simplified + D)
dotprint(transfer_func_simplified)
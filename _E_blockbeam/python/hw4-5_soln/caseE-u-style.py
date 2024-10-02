# %% [markdown]
# # Linearization: Design Study

# %%
import sympy
from sympy import *
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

init_printing()

def dotprint(expr):
    display(Math(vlatex(expr)))

# %%
t = symbols('t')
# Generalized coordinates
z, theta = symbols(r'z, \theta', cls=Function)
z = z(t)
theta = theta(t)

z_dot = z.diff(t)
theta_dot = theta.diff(t)

z_ddot = z.diff(t,2)
theta_ddot = theta.diff(t,2)

m1, m2, ell, g, F = symbols(r'm_1, m_2, \ell, g, F', real=True)

# %%
# Equations of motion in state variable form
f_of_x_and_u = Matrix([
    z_dot,
    theta_dot,
    1/m1*(-m1*g*sin(theta) + m1*z*theta_dot**2),
    1/(m2*ell**2/3 + m1*z**2)*(ell*F*cos(theta) - 2*m1*z*z_dot*theta_dot - m1*g*z*cos(theta) - m2*g*ell/2*cos(theta))
])

dotprint(f_of_x_and_u)

# %% [markdown]
# ## Deriving Nonlinear State Space Equations

# %%
state = Matrix([z, theta, z_dot, theta_dot])
dotprint(state)

# %%
state_deriv = f_of_x_and_u
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
u_eq = m1*g/ell*z_e + m2*g/2
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
A = f_of_x_and_u.jacobian(state)
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
B = f_of_x_and_u.jacobian(Matrix([F]))
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
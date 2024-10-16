# %% [markdown]
# # Ball on Beam Linearization - Individual Terms (Fixed)
# %%
import sympy as sp
from sympy import symbols, Function, diff, Matrix, sin, cos, solve, Eq
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

def dotprint(expr):
    display(Math(vlatex(expr)))

# %%
# Define symbols and functions
t = symbols('t')
z, theta = symbols('z theta', cls=Function)
z = z(t)
theta = theta(t)
z_dot = z.diff(t)
theta_dot = theta.diff(t)
z_ddot = z.diff(t, 2)
theta_ddot = theta.diff(t, 2)
m1, m2, l, g, F = symbols('m_1 m_2 l g F', real=True)

# %%
# Define nonlinear terms
nonlinear_terms = {
    'z_theta_dot_squared': z * theta_dot**2,
    'sin_theta': sin(theta),
    'z_z_dot_theta_dot': z * z_dot * theta_dot,
    'z_cos_theta': z * cos(theta),
    'cos_theta': cos(theta)
}

# %%
# Define specific equilibrium point
z_e, F_e = symbols('z_e F_e')
theta_e = 0  # We know theta_e is 0 at equilibrium
eq_point = {z: z_e, theta: theta_e, z_dot: 0, theta_dot: 0}

# %% [markdown]
# ## Linearization of Individual Terms
# %%
# Linearize each nonlinear term
linearized_terms = {}

for term_name, term in nonlinear_terms.items():
    # Compute Taylor expansion around equilibrium point
    taylor_expansion = term.series(z, z_e, 1).removeO().series(theta, theta_e, 1).removeO()
    
    # Substitute equilibrium values
    linearized = taylor_expansion.subs(eq_point)
    
    # Store the linearized term
    linearized_terms[term_name] = linearized

    print(f"Linearized {term_name}:")
    dotprint(linearized)
    print("\n")

# %% [markdown]
# ## Equilibrium Analysis
# %%
# Equations of motion
eq1 = m1 * z_ddot - m1 * nonlinear_terms['z_theta_dot_squared'] + m1 * g * nonlinear_terms['sin_theta']
eq2 = (m2*l**2/3 + m1*z**2) * theta_ddot + 2*m1*nonlinear_terms['z_z_dot_theta_dot'] + \
      m1*g*nonlinear_terms['z_cos_theta'] + m2*g*l/2*nonlinear_terms['cos_theta'] - l*F*nonlinear_terms['cos_theta']

# Set all derivatives to zero for equilibrium
eq_conditions = {z_dot: 0, theta_dot: 0, z_ddot: 0, theta_ddot: 0, theta: theta_e}
eq_eqs = [eq.subs(eq_conditions) for eq in [eq1, eq2]]

# Solve for F_e
F_e_sol = solve(eq_eqs[1], F)[0]

print("Equilibrium force F_e:")
dotprint(F_e_sol)

# %% [markdown]
# ## Linearized Equations of Motion
# %%
# Substitute linearized terms into equations of motion
eq1_lin = m1 * z_ddot - m1 * linearized_terms['z_theta_dot_squared'] + m1 * g * linearized_terms['sin_theta']
eq2_lin = (m2*l**2/3 + m1*z_e**2) * theta_ddot + 2*m1*linearized_terms['z_z_dot_theta_dot'] + \
          m1*g*linearized_terms['z_cos_theta'] + m2*g*l/2*linearized_terms['cos_theta'] - l*F*linearized_terms['cos_theta']

print("Linearized equation 1:")
dotprint(eq1_lin)
print("\nLinearized equation 2:")
dotprint(eq2_lin)

# %% [markdown]
# ## State Space Representation
# %%
# Define state vector
x = Matrix([z - z_e, theta - theta_e, z_dot, theta_dot])

# Define input
u = Matrix([F - F_e_sol])

# Solve for z_ddot and theta_ddot
solved = solve([eq1_lin, eq2_lin], (z_ddot, theta_ddot))

# Define state equations
dx_dt = Matrix([
    z_dot,
    theta_dot,
    solved[z_ddot],
    solved[theta_ddot]
])

# Compute A and B matrices
A = dx_dt.jacobian(x)
B = dx_dt.jacobian(u)

print("A matrix:")
dotprint(A)
print("\nB matrix:")
dotprint(B)

# %% [markdown]
# ## Transfer Function
# %%
C = Matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
D = Matrix([[0], [0]])
s = symbols('s')

transfer_func = C * (s*sp.eye(4) - A).inv() * B + D
transfer_func_simplified = sp.simplify(transfer_func)

print("Transfer function:")
dotprint(transfer_func_simplified)
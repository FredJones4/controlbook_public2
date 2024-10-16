# %% [markdown]
# # Mass-Spring-Damper Linearization

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
z = symbols('z', cls=Function)
z = z(t)
z_dot = z.diff(t)
z_ddot = z.diff(t,2)

m, b, k, F = symbols('m, b, k, F', real=True)

# %%
LHS = m * z_ddot + b * z_dot + k * z
RHS = F
dynamics = Eq(LHS, RHS)
dotprint(dynamics)

# %% [markdown]
# ## Deriving Nonlinear State Space Equations
# 
# Solve for our highest derivative, $\ddot{z}$:

# %%
solve_dict = solve(dynamics, (z_ddot,), simplify=True, dict=True)[0]
dotprint(solve_dict)

# %% [markdown]
# Create our state vector $x = [x_1, x_2]$:

# %%
state = MatrixSymbol('x', 2, 1)
dotprint(state)

# %% [markdown]
# We have $x = [x_1, x_2] = [z, \dot{z}]$, so $\dot{x} = [\dot{z}, \ddot{z}]$. Let's put that in a vector:

# %%
z_ddot_expr = solve_dict[z_ddot]
state_deriv = Matrix([z_dot, z_ddot_expr])
dotprint(state_deriv)

# %% [markdown]
# Now we can substitute in our $x_1, x_2$ values:

# %%
# Dictionary for substitutions
subs_dict = {
    z: state[0],
    z_dot: state[1],
}
f_expr = state_deriv.subs(subs_dict)
dotprint(f_expr)

# %% [markdown]
# We now have our state space equations. Note that these are already linear, so we don't need to perform linearization.
# 
# ## Equilibrium Points
# 
# Let's find the equilibrium points:

# %%
# This equation represents f(x,u) = 0
equilibrium_equation = Eq(f_expr, Matrix([0,0]))
eq_solve_dict = solve(equilibrium_equation, (state[0], state[1], F), simplify=True, dict=True)[0]
dotprint(eq_solve_dict)

# %% [markdown]
# From this, we can see that at equilibrium:
# 
# 1. $z_e$ can be any value
# 2. $\dot{z}_e = 0$
# 3. $F_e = kz_e$
# 
# The system is already linear!

# %% [markdown]
# ## A and B Matrices
# 
# Even though the system is already linear, let's derive the A and B matrices for completeness:

# %%
A = f_expr.jacobian(state)
dotprint(A)

# %%
B = f_expr.jacobian([F])
dotprint(B)

# %% [markdown]
# Note that A and B are constant matrices and do not depend on the state or input. This is characteristic of linear systems.

# %% [markdown]
# ## Conclusion
# 
# The equations of motion for this system are linear and do not require linearization. However, the techniques of feedback linearization can still be applied if we want to achieve a specific equilibrium output.
# 
# For the mass-spring-damper system, a non-zero force is required to hold the mass at any non-zero equilibrium position (i.e., $F = kz + F_0$), due to the spring force generated when the position of the mass is not zero.
# %%
#Import Relevant Packages
import sympy
from sympy import *

#Initialize LaTeX printing
sympy.init_printing()

from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

init_printing()

#Define the dotprint to get '.' above variables
def dotprint(expr):
    display(Math(vlatex(expr)))

# %%
#Define configuration variables and their derivatives
t = symbols('t')

theta = symbols(r'\theta', cls=Function)
theta = theta(t)

z = symbols(r'z',cls=Function)
z = z(t)

theta_dot = theta.diff(t)
z_dot = z.diff(t)

theta_ddot = theta.diff(t,2)
z_ddot = z.diff(t,2)

dotprint(z_dot),dotprint(theta_dot)

# %%
#Define all other symbols
m1, m2, ell, Po, g, F, b, tau = symbols(r'm_1, m_2, \ell, P_o, g, F, b \tau', real=True)

# %%
#Define position matrix
P1 = Matrix([z*cos(theta),z*sin(theta), 0 ])
P2 = Matrix([(ell/2)*cos(theta), (ell/2)*sin(theta), 0])
P1,P2

# %%
#Differentiate to get velocity
V1 = P1.diff(t)
V2 = P2.diff(t)
dotprint(V1),dotprint(V2)

# %%
#Derive linear kinetic energy
K_lin = Rational(1,2)*m1*V1.T @ V1 + Rational(1,2)*m2*V2.T @ V2

dotprint(K_lin)

# %%
#Derive Rotational Kinetic Energy
w = Matrix([0,0,theta_dot])
J = sympy.zeros(3,3)
J[2,2] = m2*ell**2/12
K_rot = Rational(1,2)*w.T @ J @ w
dotprint(K_rot)

# %%
#Combine for total Kinetic Energy
K_total = K_lin + K_rot
K = simplify(K_total)
dotprint(K)

# %%
#Derive Potential Energy
P = Po + Rational(1,2)*(m2)*g*ell*sin(theta) + m1*g*z*sin(theta)
dotprint(P)

# %%
#Define Lagrangian
L = K[0] - P
dotprint(L)

# %%
#Define input forces tau
tautemp = Matrix([[0],[F*ell*cos(theta)]])
RHS1 = Matrix([tau])
LHS1 = Matrix([tautemp]) # changed tautemp for tau : NOTE: this might need
dynamics = Eq(LHS1,RHS1,evaluate=False)

dotprint(dynamics)
LHS1 = RHS1

# %%
#Non-Conservative Forces
B = Matrix([[0,0],[0,0]]) #What is b in this scenario?
q = Matrix([z, theta])
q_dot = Matrix([z_dot,theta_dot])
Bq_dot = B @ q_dot
F_total = tautemp - Bq_dot

dotprint(F_total)

# %%
#Derive Langrangian function
def derive_Lagrangian(L, q, q_dot):
    term_1 = (sympy.tensor.derive_by_array(L, q_dot)).diff(t)
    term_2 = sympy.tensor.derive_by_array(L, q)
    return term_1 - term_2

# %%
#Assemble Final Equation
LHS = Matrix(derive_Lagrangian(L,q,q_dot))
RHS = Matrix(tautemp - Bq_dot)
Euler_Lagrange = Eq(LHS, RHS, evaluate=False)

dotprint(Euler_Lagrange)

# %%
solve_dict = solve(Euler_Lagrange,[z_ddot,theta_ddot], simplify = True, dict = True)[0]
dotprint(solve_dict)

# %%
state = MatrixSymbol('x', 4,1)
dotprint(state)

# %%
theta_ddot_expr = solve_dict[theta_ddot]
z_ddot_expr = solve_dict[z_ddot]
state_deriv = Matrix([theta_dot, z_dot, theta_ddot_expr, z_ddot_expr])
dotprint(state_deriv)

# %%
dotprint(theta_ddot_expr)

# %%
#Dictionary for substitutions
subs_dict = {
    theta: state[0],
    z: state[1],
    theta_dot: state[2],
    z_dot: state[3],
}

f_expr = state_deriv.subs(subs_dict)
dotprint(f_expr)

# %%
#Equilibrum equation and solving

equilibrium = Eq(f_expr, Matrix([0,0,0,0]))
eq_solve_dict = solve(equilibrium, (state[0],state[1],state[2],state[3], tau), simplify = True, dict = True)[0]
dotprint(eq_solve_dict)

# %%
A = f_expr.jacobian(state)
dotprint(A)

# %%
dotprint(f_expr)

# %%
force_term = F*ell*cos(theta)
B = f_expr.jacobian([force_term])
dotprint((B))

# %%

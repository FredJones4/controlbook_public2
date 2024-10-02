# %%
import sympy
from sympy import *
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

init_printing()

def dotprint(expr):
    display(Math(vlatex(expr)))

# %%
y = Matrix([[1,1]])

dotprint(y[0])

# %%
t = symbols('t')
z, h, theta = symbols(r'z, h, \theta', cls=Function)
z = z(t)
theta = theta(t)
h = h(t)
mc, mm, jc, d, Fl, Fr, g, mu= symbols(r'mc, mm, jc, d, Fl, Fr, g, \mu')
z_dot = z.diff(t)
theta_dot = theta.diff(t)
h_dot = h.diff(t)

z_ddot = Derivative(z, (t,2))
theta_ddot = Derivative(theta, (t,2))
h_ddot = Derivative(h, (t,2))

dotprint(mu)


# %%
q = Matrix([z,h, theta])
q_dot = q.diff(t)

pc = Matrix([z,h])
pl = Matrix([z-d*cos(theta), h-d*sin(theta)])
pr = Matrix([z+d*cos(theta), h+d*sin(theta)])

vc = pc.diff(t)
vl = pl.diff(t)
vr = pr.diff(t)

Kc = Rational(1,2)*mc*vc.T @ vc +Matrix([Rational(1,2)*theta_dot**2*jc])
Kl = Rational(1,2)*mm*vl.T @ vl
Kr = Rational(1,2)*mm*vr.T @ vr


K = Kc+Kl+Kr


dotprint(K[0])


# %%
Pc = mc*g*pc[1]
Pr = mm*g*pr[1]
Pl = mm*g*pl[1]

P = Pc+Pr+Pl

dotprint(P)

# %%
L = K[0]-P

dotprint(L)



# %%
tau = Matrix([[-(Fl+Fr)*sin(theta)],[(Fl+Fr)*cos(theta)], [(Fl-Fr)*d]])
b = Matrix([[mu,0,0],[0,0,0],[0,0,0]])

bqdot = b@q_dot

F_total = tau +bqdot

dotprint(F_total)

# %%
# Modified derive_Lagrangian function
def derive_Lagrangian(L, q, q_dot):
    # Term 1: Derivative of d(L)/d(q_dot)
    dL_dqdot = Matrix([L.diff(qi_dot) for qi_dot in q_dot])
    term_1 = dL_dqdot.diff(t)
    
    # Term 2: Derivative of d(L)/d(q)
    dL_dq = Matrix([L.diff(qi) for qi in q])
    term_2 = dL_dq
    
    # Return the Euler-Lagrange equation
    return term_1 - term_2


# %%
LHS = derive_Lagrangian(L, q, q_dot)

RHS = F_total

Euler_Lagrange = Eq(LHS, RHS)

dotprint(simplify(Euler_Lagrange))

desired = Matrix([[z_ddot],[h_ddot], [theta_ddot]])
solution = solve(Euler_Lagrange, desired)

dotprint(desired)

print("Solution is")
dotprint(solution)

#dotprint(solve(Euler_Lagrange, desired))

#



# %%

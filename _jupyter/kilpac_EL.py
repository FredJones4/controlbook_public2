#%%
import sympy as sp
from IPython.display import Math, display

# if displaying is not working, you can try uncommenting the following lines: 
# from sympy import init_printing
# init_printing(use_latex="mathjax", latex_mode="equation")


#%%
####################################################################################
#example from Case Study A
####################################################################################
# defining all necessary mathematical variables as sympy "symbols" so that we can use them accordingly
t, m, ell, g = sp.symbols('t, m, ell, g')

# define theta as a function of time so we can take the time derivative
theta = sp.Function('theta')(t)

#%%
#calculate the kinetic energy
K = 1/2.0*m*ell**2/3*theta.diff(t)**2

#calculate the potential energy
P = m*g*ell/2.0*sp.sin(theta)

#calculate the lagrangian
L = K-P

#%%
#calculate the Euler-Lagrange equations of motion for theta variable
EL_case_studyA = sp.diff(sp.diff(L,sp.diff(theta,t)), t) - sp.diff(L,theta)

# previous lines to calculate EL_case_studyA could also be accomplished like this:
# theta_dot = theta.diff(t)
# EL_case_studyA = (L.diff(theta_dot)).diff(t) - L.diff(theta)


print("\n\n\n\n Euler-Lagrange equations for case study A:")
sp.pretty_print(EL_case_studyA)

display(EL_case_studyA)


#%%
####################################################################################
#example from Case Study B (with vectors/matrices)
####################################################################################
# importing these functions directly to make life a little easier and code a little more readable
from sympy import sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, latex, init_printing
#init_printing(use_latex="mathjax", latex_mode="equation")

#defining mathematical variables (called symbols in sp) and time varying functions like z and theta
t, m1, m2, ell, g = symbols('t, m1, m2, ell, g')
z = Function('z')(t)
theta = Function('theta')(t)

#defining generalized coords and derivatives
q = Matrix([[z], [theta]])
qdot = diff(q, t)  #q.diff(t)

#defining the kinetic energy
p1 = Matrix([[z+ell/2.0*sin(theta)], [ell/2.0*cos(theta)], [0]])
p2 = Matrix([[z], [0], [0]])

v1 = diff(p1, t)
v2 = diff(p2, t)

omega = Matrix([[0], [0], [diff(theta,t)]])
J = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, m1*ell**2/12.0]])

#%%
K = simplify(0.5*m1*v1.T*v1 + 0.5*m2*v2.T*v2 + 0.5*omega.T*J*omega)
K = K[0,0]

print("\n\n\n\n Kinetic energy for case study B (compare with book)")
display(K)

#%%
#defining potential energy (MUST BE A MATRIX as well to do L = K-P)
P = m1*g*ell/2.0*(cos(theta)-1)

#calculate the lagrangian, using simplify intermittently can help the equations to be
#simpler, there are also options for factoring and grouping if you look at the sympy
#documentation.
L = simplify(K-P)

#%%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like -B*q_dot and tau)
EL_case_studyB = simplify( diff(diff(L, qdot), t) - diff(L, q) )

# if you want to use the latex output, you can put it in a latex editor, or website like
# this - https://arachnoid.com/latex/. Some interactive python IDEs will render the latex
# code directly (such as Jupyter notebooks), but I can't get pycharm to do it.
print("\n\n\n\n Latex code:\n", latex(EL_case_studyB))

# we can write this latex code to file if that's helpful:
out_file = open("latex_EL_case_studyB.txt","w")
out_file.write(latex(EL_case_studyB))
out_file.close()

# we can also just see the output like this:
print("\n\n\n\n Solution for case study B (compare with book solution)")
display(EL_case_studyB)


# %%
####################################################################################
#re-do Case Study A using "dynamicsymbols" that are by default functions of time, and
# printing those variables a little differently (with dot notation for derivatives)
####################################################################################
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex

# defining all necessary mathematical variables as sympy "symbols" so that we can use them accordingly
t, m, ell, g = sp.symbols('t, m, ell, g')

# define theta and theta_dot as functions of time so we can take the time derivative
theta = dynamicsymbols('theta')
theta_dot = theta.diff(t)

#calculate the kinetic energy
K = 1/2.0*m*ell**2/3*theta_dot**2

#calculate the potential energy
P = m*g*ell/2.0*sp.sin(theta)

#calculate the lagrangian
L = K-P

#calculate the Euler-Lagrange equations of motion for theta variable
EL_case_studyA_v2 = sp.simplify((L.diff(theta_dot)).diff(t) - L.diff(theta))

print("\n\n\n\n Euler-Lagrange equations for case study A:")
display(EL_case_studyA_v2)

# to see the dot notation:
display(Math(vlatex(EL_case_studyA_v2)))

# if you want dot notation (instead of d/dt) in your latex code, this will work too:
vlatex(EL_case_studyA_v2)


# %%

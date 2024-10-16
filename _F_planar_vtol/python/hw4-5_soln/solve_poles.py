import sympy as sp
import control
import matplotlib.pyplot as plt
import numpy as np

# Define symbols
s, mc, mr, Jc, d, Fe, mu = sp.symbols('s mc mr Jc d Fe mu')
F, tau, H, Z, Theta = sp.symbols('F tau H Z Theta')

def longitudinal_dynamics():
    # Equation (1) in s-domain
    eq = sp.Eq((mc + 2*mr) * s**2 * H, F)
    
    # Solve for H(s)
    H_s = sp.solve(eq, H)[0]
    
    # Transfer function H(s)/F(s)
    TF_H_F = sp.simplify(H_s / F)
    
    return TF_H_F

def lateral_dynamics():
    # Define a1 and a2
    a1 = mc + 2*mr
    a2 = Jc + 2*mr*d**2
    
    # Equations (3) and (4) in s-domain
    eq1 = sp.Eq(a1 * s**2 * Z + mu * s * Z, -Fe * Theta)
    eq2 = sp.Eq(a2 * s**2 * Theta, tau)
    
    # Solve for Z and Theta in terms of tau
    solution = sp.solve((eq1, eq2), (Z, Theta))
    
    # Transfer functions
    TF_Z_tau = sp.simplify(solution[Z] / tau)
    TF_Theta_tau = sp.simplify(solution[Theta] / tau)
    
    # Transfer function from Theta to Z
    TF_Z_Theta = sp.simplify(-Fe / (a1 * s**2 + mu * s))
    
    return TF_Z_tau, TF_Theta_tau, TF_Z_Theta

def print_results(TF_H_F, TF_Z_tau, TF_Theta_tau, TF_Z_Theta):
    print("Longitudinal Dynamics:")
    print("H(s)/F(s) = ")
    sp.pprint(TF_H_F)
    
    print("\nLateral Dynamics:")
    print("Z(s)/tau(s) = ")
    # NOTE:
    """
    The answer is equivalent to the online PDF, but the code output is more expanded.
    In the PDF, a1 = mc + 2mr and a2 = Jc + 2mr*d^2,
    which when expanded matches the code output.
    """
    sp.pprint(TF_Z_tau)
    print("\nTheta(s)/tau(s) = ")
    """
    PDF: Theta(s)/tau(s) = (1 / a2) * (1 / s^2)
    Code: Theta(s)/tau(s) = 1 / (s^2 * (Jc + 2d^2mr))
    These are equivalent.
    """
    sp.pprint(TF_Theta_tau)
    print("\nZ(s)/Theta(s) = ")
    """
    PDF: Z(s)/Theta(s) = (-F0 / a1) * (1 / (s * (s + (μ / a1))))
    Code: Z(s)/Theta(s) = -Fe / (s*(μ + s*(mc + 2*mr)))
    These are equivalent. The code output is correct, with a1 expanded to (mc + 2mr).
    """
    sp.pprint(TF_Z_Theta)

def draw_block_diagrams():
    # This function would typically use a library like networkx or graphviz
    # to draw the block diagrams. For simplicity, we'll just print text representations.
    print("\nLongitudinal Dynamics Block Diagram (text representation):")
    print("F(s) -> [1 / ((mc + 2mr) * s^2)] -> H(s)")
    
    print("\nLateral Dynamics Block Diagram (text representation):")
    print("tau(s) -> [1 / (Jc + 2mr*d^2) * s^2] -> Theta(s) -> [-Fe / ((mc + 2mr)*s^2 + mu*s)] -> Z(s)")

if __name__ == "__main__":
    # (a) and (b) Longitudinal dynamics
    TF_H_F = longitudinal_dynamics()
    
    # (c) Lateral dynamics
    TF_Z_tau, TF_Theta_tau, TF_Z_Theta = lateral_dynamics()
    
    # Print results
    print_results(TF_H_F, TF_Z_tau, TF_Theta_tau, TF_Z_Theta)
    
    # (d) Block diagrams
    draw_block_diagrams()
import sympy as sp
import control
import matplotlib.pyplot as plt
import numpy as np

# Define symbols
s, m1, m2, l, g, ze, A = sp.symbols('s m1 m2 l g ze A')
F, Z, Theta = sp.symbols('F Z Theta')

def equations_of_motion():
    # Linearized equations of motion
    eq1 = sp.Eq(m1 * s**2 * Z + m1 * g * Theta, 0)
    eq2 = sp.Eq(m1 * g * Z + A * s**2 * Theta, l * F)
    
    # Define A
    A_def = sp.Eq(A, (m2 * l**2 / 3) + m1 * ze**2)
    
    return eq1, eq2, A_def

def transfer_functions():
    eq1, eq2, A_def = equations_of_motion()
    
    # Solve for Z and Theta in terms of F
    solution = sp.solve((eq1, eq2), (Z, Theta))
    
    # Transfer functions
    TF_Z_F = sp.simplify(solution[Z] / F)
    TF_Theta_F = sp.simplify(solution[Theta] / F)
    
    # Transfer function from Theta to Z
    TF_Z_Theta = sp.simplify(-g / s**2)
    
    return TF_Z_F, TF_Theta_F, TF_Z_Theta

def simplified_transfer_functions():
    # Neglect m1*g*Z term in eq2
    eq1 = sp.Eq(m1 * s**2 * Z + m1 * g * Theta, 0)
    eq2_simplified = sp.Eq(A * s**2 * Theta, l * F)
    
    # Solve for Theta in terms of F
    Theta_F = sp.solve(eq2_simplified, Theta)[0] / F
    
    # Use the relationship from eq1 to find Z in terms of Theta
    Z_Theta = sp.solve(eq1, Z)[0] / Theta
    
    # Combine to get Z in terms of F
    Z_F = sp.simplify(Z_Theta * Theta_F)
    
    return Z_F, Theta_F, Z_Theta

def print_results(TF_Z_F, TF_Theta_F, TF_Z_Theta, simplified=False):
    print("Transfer Functions:")
    print(f"{'Simplified ' if simplified else ''}Z(s)/F(s) = ")
    sp.pprint(TF_Z_F)
    print(f"\n{'Simplified ' if simplified else ''}Theta(s)/F(s) = ")
    sp.pprint(TF_Theta_F)
    print(f"\n{'Simplified ' if simplified else ''}Z(s)/Theta(s) = ")
    sp.pprint(TF_Z_Theta)

def draw_block_diagram():
    # This function would typically use a library like networkx or graphviz
    # to draw the block diagram. For simplicity, we'll just print a text representation.
    print("\nBlock Diagram (text representation):")
    print("F(s) -> [l/(A*s^2)] -> Theta(s) -> [-g/s^2] -> Z(s)")

if __name__ == "__main__":
    # (a) and (b) Full transfer functions
    TF_Z_F, TF_Theta_F, TF_Z_Theta = transfer_functions()
    print_results(TF_Z_F, TF_Theta_F, TF_Z_Theta)
    
    # (c) Simplified transfer functions
    TF_Z_F_simple, TF_Theta_F_simple, TF_Z_Theta_simple = simplified_transfer_functions()
    print("\nSimplified Transfer Functions:")
    print_results(TF_Z_F_simple, TF_Theta_F_simple, TF_Z_Theta_simple, simplified=True)
    
    # (d) Block diagram
    draw_block_diagram()
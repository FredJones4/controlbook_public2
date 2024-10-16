import sympy as sp

# Define symbolic variables
s, m, b, k, kP, kD = sp.symbols('s m b k kP kD')

def explain_and_print(expr, explanation):
    print(f"{explanation}:")
    print(f"{expr}\n")

# Part (a): Open loop transfer function and poles
def open_loop_analysis():
    # Define the open loop transfer function
    G_s = 1 / (m * s**2 + b * s + k)
    explain_and_print(G_s, "Open loop transfer function")

    # Substitute the given values
    G_s_substituted = G_s.subs({m: 5, b: 0.5, k: 3})
    explain_and_print(G_s_substituted, "Open loop transfer function with substituted values")

    # Find the open loop poles
    denominator = sp.denom(G_s_substituted)
    open_loop_poles = sp.solve(denominator, s)
    explain_and_print(open_loop_poles, "Open loop poles")

# Part (b): Closed loop transfer function and poles
def closed_loop_analysis():
    # Define the closed loop transfer function
    G_s = 0.2 / (s**2 + 0.1*s + 0.6)
    H_s = kP + kD * s
    T_s = (G_s * H_s) / (1 + G_s * H_s)
    explain_and_print(T_s, "Closed loop transfer function")

    # Simplify the closed loop transfer function
    T_s_simplified = sp.simplify(T_s)
    explain_and_print(T_s_simplified, "Simplified closed loop transfer function")

    # Extract the characteristic equation
    char_eq = sp.denom(T_s_simplified)
    explain_and_print(char_eq, "Characteristic equation")

# Part (c): Find kP and kD for desired poles
def find_controller_gains():
    # Define the desired characteristic equation
    desired_char_eq = s**2 + 2.5*s + 1.5

    # Equate coefficients
    eq1 = sp.Eq(0.1 + 0.2*kD, 2.5)
    eq2 = sp.Eq(0.6 + 0.2*kP, 1.5)

    # Solve for kP and kD
    solution = sp.solve((eq1, eq2), (kP, kD))
    explain_and_print(solution, "Controller gains (kP, kD)")

# Execute the analysis
if __name__ == "__main__":
    print("Part (a): Open Loop Analysis")
    open_loop_analysis()

    print("Part (b): Closed Loop Analysis")
    closed_loop_analysis()

    print("Part (c): Finding Controller Gains")
    find_controller_gains()
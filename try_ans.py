import sympy as sp
import control
import matplotlib.pyplot as plt
import numpy as np

# Define symbols
s, m, b, k = sp.symbols('s m b k')
kp, kd = sp.symbols('kp kd')

# (a) Open loop transfer function and poles
def open_loop_tf_and_poles():
    # Transfer function
    G = 1 / (m * s**2 + b * s + k)
    
    # Substitute values
    G_subs = G.subs({m: 5, b: 0.5, k: 3})
    
    # Convert to fraction
    G_frac = sp.fraction(G_subs)
    
    # Find poles
    poles = sp.solve(G_frac[1], s)
    
    return G_subs, poles

# (b) Closed loop transfer function and poles
def closed_loop_tf_and_poles():
    G = 0.2 / (s**2 + 0.1*s + 0.6)
    H = kp + kd*s
    
    # Closed loop transfer function
    CL = G * H / (1 + G * H)
    CL_simplified = sp.simplify(CL)
    
    # Convert to fraction
    CL_frac = sp.fraction(CL_simplified)
    
    # Closed loop poles (roots of the denominator)
    poles = sp.solve(CL_frac[1], s)
    
    return CL_simplified, poles

# (c) Solve for kp and kd
def solve_for_gains():
    # Desired characteristic equation: s^2 + 2.5s + 1.5
    eq1 = sp.Eq(0.1 + 0.2*kd, 2.5)
    eq2 = sp.Eq(0.6 + 0.2*kp, 1.5)
    
    solution = sp.solve((eq1, eq2), (kp, kd))
    return solution

# (d) Simulate step response
def simulate_step_response(kp, kd):
    num = [0.2 * kp]
    den = [1, 0.1 + 0.2*kd, 0.6 + 0.2*kp]
    sys = control.TransferFunction(num, den)
    
    t = np.linspace(0, 10, 1000)
    t, y = control.step_response(sys, T=t)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t, y)
    plt.title('Step Response')
    plt.xlabel('Time (s)')
    plt.ylabel('Position (m)')
    plt.grid(True)
    plt.show()

# Main execution
if __name__ == "__main__":
    # (a) Open loop analysis
    G_ol, ol_poles = open_loop_tf_and_poles()
    print("Open Loop Transfer Function:")
    sp.pprint(G_ol)
    print("Open Loop Poles:")
    sp.pprint(ol_poles)
    
    # (b) Closed loop analysis
    CL_tf, cl_poles = closed_loop_tf_and_poles()
    print("\nClosed Loop Transfer Function:")
    sp.pprint(CL_tf)
    print("Closed Loop Poles (as functions of kp and kd):")
    sp.pprint(cl_poles)
    
    # (c) Solve for kp and kd
    gains = solve_for_gains()
    print("\nSolved gains:")
    print(f"kp = {gains[kp]}")
    print(f"kd = {gains[kd]}")
    
    # (d) Simulate step response
    kp_value = float(gains[kp])
    kd_value = float(gains[kd])
    simulate_step_response(kp_value, kd_value)
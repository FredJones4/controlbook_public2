#%%
# NOTE: This can be run as a Jupyter Notebook in VSCode.
import sympy as sp
def calculate_system_type_by_denominator(den_expr: sp.core.mul.Mul):
    """
    Optimized system type calculation that directly checks for multiples of s in the denominator.
    """
    # Define the Laplace variable 's'
    s = sp.symbols('s')
    # print("tf:", den_expr)
    simplified_den = den_expr
    # Step 3: Use sp.roots() to get the poles and their multiplicities
    poles_with_multiplicities = sp.roots(simplified_den, s)
    # Step 4: Convert the dictionary of poles and multiplicities into a list that includes all multiples
    poles = []
    for pole, multiplicity in poles_with_multiplicities.items():
        poles.extend([pole] * multiplicity)  # Add the pole multiple times based on its multiplicity
        # Count the number of integrators (poles at the origin, s = 0)
    num_integrators = sum(1 for pole in poles if pole.is_zero)

    # Determine the system type based on the number of integrators
    system_type = num_integrators

    # Define the steady-state error expectations based on the system type
    if system_type == 0:
        step_error = "Finite error (calculate)"
        ramp_error = "Infinite error"
        parabola_error = "Infinite error"
    elif system_type == 1:
        step_error = "Zero error"
        ramp_error = "Finite error (calculate)"
        parabola_error = "Infinite error"
    elif system_type == 2:
        step_error = "Zero error"
        ramp_error = "Zero error"
        parabola_error = "Finite error (calculate)"
    else:
        step_error = "Zero error"
        ramp_error = "Zero error"
        parabola_error = "Zero error"

    # Provide recommendations based on system type
    recommendations = {
        0: "Consider adding integrators to improve tracking for ramps and parabolas.",
        1: "System tracks steps perfectly. Consider adding another integrator for ramp tracking.",
        2: "System tracks steps and ramps perfectly. Add more integrators if parabola tracking is needed.",
        3: "System tracks steps, ramps, and parabolas perfectly."
    }

    # Output the results in a dictionary format
    result = {
        "System Type": system_type,
        "Step Input Error": step_error,
        "Ramp Input Error": ramp_error,
        "Parabola Input Error": parabola_error,
        "Recommendation": recommendations.get(system_type, "No further improvements needed.")
    }

    return result
def tracking_error(tf: sp.core.mul.Mul, q: int):
    """
    Calculate the steady-state tracking error for a system based on the transfer function
    and the input type q (step, ramp, parabola, etc.).
    
    Args:
    tf (sympy expression): The transfer function as a sympy expression.
    q (int): The maximum input type. For example, 0 for step, 1 for ramp, 2 for parabola, etc.
    
    Returns:
    dict: A dictionary containing the steady-state tracking error for each q_cur from 0 to q.
    """
    # Define the Laplace variable 's'
    s = sp.symbols('s')
    s, k_d, k_p, k_i, ell, m, b, m1, m2, g, A, Js, Jp, k = sp.symbols('s k_d k_p k_i ell m b m1 m2 g A Js Jp k')
    # Initialize an empty dictionary to store the tracking errors for each q_cur
    tracking_errors = {}

    # Iterate over all values of q_cur from 0 to q (inclusive)
    for q_cur in range(q + 1):
        # Tracking error formula: 1 / ((s**q_cur) * tf)
        tracking_error_expr = 1 / ((s**q_cur) * (1+ tf)) #TODO: check for accuracy
        
        # Take the limit as s approaches 0 (final steady-state error)
        final_tracking_error = sp.limit(tracking_error_expr, s, 0)
        
        # Store the result in the dictionary with q_cur as the key
        tracking_errors[q_cur] = final_tracking_error
    
    return tracking_errors
def steady_state_error(tf: sp.core.mul.Mul, d: sp.core.mul.Mul):
    s = sp.symbols('s')
    return sp.simplify(sp.limit((tf)*d,s,0))
def find_q_for_steady_state_error(P,C):
    """
    This function takes a transfer function `tf` and iteratively increases `q` to compute the steady-state error.
    It returns the smallest `q` where the steady-state error is finite and non-zero, along with the steady-state error.
    
    Args:
    tf (sympy expression): The closed-loop transfer function.
    A (sympy symbol or value): The amplitude of the disturbance.
    
    Returns:
    tuple: (int, sympy expression) The smallest value of q and the corresponding steady-state error, 
           where the steady-state error is finite and non-zero.
    """
    # Define necessary symbols
    A = sp.symbols('A')
    s = sp.symbols('s')
    q = sp.symbols('q', positive=True, integer=True)
    
    # Start with q = 0 and increment
    for q_val in range(10):  # Arbitrary upper limit of 10, you can adjust this if needed
        
        # Step 4: Define disturbance D(s) = A / s^(q+1)
        D = A / s**(q + 1)
        
        # Step 5: Calculate the error transfer function E(s) = (1 / (1 + tf)) * D
        error_tf = s*(P / (1 + P*C)) * D
        # Step 6: Substitute q = q_val before taking the limit
        error_tf_q = error_tf.subs(q, q_val)
        
        # Step 7: Simplify the expression
        simplified_error_tf = sp.simplify(error_tf_q)
        
        # Step 8: Compute the steady-state error by taking the limit as s -> 0
        steady_state_error = sp.limit(simplified_error_tf, s, 0)
        
        # Step 9: Check if the steady-state error is valid (not zero, not infinity)
        if steady_state_error != 0 and not steady_state_error.has(sp.oo, -sp.oo):
            return q_val, steady_state_error
    
    # If no valid q is found, return None (or raise an exception depending on preference)
    return None, None
# %%
# from sympy_system_type import *
def analyze_transfer_function(P):
    # Define the symbols used in the analysis
    s = sp.symbols('s')
    k_d, k_p, k_i,  = sp.symbols('k_d k_p k_i')

    # PD and PID controllers
    C_pid = (k_d*s**2 + k_p*s + k_i)/s
    C_pd = k_d*s + k_p

    # PD Form: System Type
    print("\n\n\033[1mQuestion A-1: When the controller for the single link robot arm is PD control, what is the system type? Characterize the steady state error when the reference input is a step, a ramp, and a parabola.\033[0m")
    transfer_function_pd = P * C_pd
    den_expr = sp.denom(transfer_function_pd)
    system_info_pd = calculate_system_type_by_denominator(den_expr)
    print(system_info_pd)

    # PD Form: Tracking Error
    q_pd = system_info_pd['System Type']
    tracking_error_pd = tracking_error(transfer_function_pd, q_pd)
    print("Tracking Error for PD form: ", tracking_error_pd)

    # PID Form: System Type
    print("\n\n\033[1mQuestion A-2: How does this change if you add an integrator?\033[0m")
    transfer_function_pid = P * C_pid
    # transfer_function_pid = sp.simplify(transfer_function_pid)
    den_expr = sp.denom(transfer_function_pid)
    system_info_pid = calculate_system_type_by_denominator(den_expr)
    print(system_info_pid)

    # PID Form: Tracking Error
    q_pid = system_info_pid['System Type']
    tracking_error_pid = tracking_error(transfer_function_pid, q_pid)
    print("Tracking error for PID form:", tracking_error_pid)

    # Question B: Steady State Error to Constant Input Disturbance
    print("\n\n\033[1mQuestion B: Consider the case where a constant disturbance acts at the input to the plant (for example gravity in this case). What is the steady state error to a constant input disturbance when the integrator is not present and when it is present?\033[0m")
    # PD Form
    tf_d = s * P / (1 + P * C_pd)
    tf_d = sp.simplify(tf_d) ## KEY FINDING, issue with the extra term otherwise
    den_expr = sp.denom(tf_d)
    q_d_pd = calculate_system_type_by_denominator(den_expr)['System Type']
    print("System type:", q_d_pd)
    # SSE - PD Form
    A = sp.symbols('A')
    D = A / (s**(q_d_pd + 1))
    sse_pd = steady_state_error(tf_d, D)
    print("PD Form, Steady State error: ", sse_pd)

    # PID Form 
    # NOTE: Out of necessity and to show a different way of handling it, a second SSE function is written.
    q_d_pid, sse_pid = find_q_for_steady_state_error(P,C_pid)
    print(f"PID q = {q_d_pid}")
    print("PID Form, Steady State error: ", sse_pid)
# %% [markdown]

# Homework Problems
#%%
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
from IPython.display import Math, display
def special_display(x):
    return display(Math(vlatex(x)))
# %% [markdown]

# Examples A-C
# %%
print("Case Study A \n\n\n\n")
s, ell, m, b = sp.symbols('s ell m b')
P = (3/(m*ell**2)) / (s*(s+3*b/(m*ell**2)))

print("Plant:")
special_display(P)
analyze_transfer_function(P)
#%%
print("Case Study B\n\n\n\n")
m1, m2, g = sp.symbols('m1 m2 g')
numerator = -1 / (m1 * (ell/6) + m2 * (2*ell/3))
denominator = s**2 - ((m1 + m2) * g) / (m1 * (ell/6) + m2 * (2*ell/3))
P_B_a =  numerator / denominator
print("Inner Loop\n\nPlant:\n")
special_display(P_B_a)
analyze_transfer_function(P_B_a)

P_B_b = (-(2*ell/3*s**2 - g)/s**2)
print("Outer Loop\n\nPlant:\n")
special_display(P_B_b)
analyze_transfer_function(P_B_b)
#%%
print("Case Study C \n\n\n\n")
print("Inner Case")
Js, Jp = sp.symbols("Js Jp")
P_prob_C_inner = 1/((Js + Jp)*s**2)
print("Plant:")
special_display(P_prob_C_inner)
analyze_transfer_function(P_prob_C_inner)
# NOTE: the output for tracking error for pd form on this outer loop is 1/kp, when the book 
# says 1/(1+kp). A more in-depth look might be necessary.
print("Outer Loop\n\n\n\n")
k = sp.symbols('k')
P_prob_C_outer = (b*s/Jp + k/Jp)/(s**2 + b/Jp*s + k/Jp)
print("Plant:")
special_display(P_prob_C_outer)
analyze_transfer_function(P_prob_C_outer)
# %% [markdown]
# Examples D, E, F are left to the reader

#%%
# print("Case Study D:\n\n\n\n")


# print("Case Study E:\n\n\n\n")

# print("Case Study F:\n\n\n\n")



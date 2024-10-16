import sympy as sp

def calculate_system_type_by_array(num, den):
    """
    Calculate the system type based on the number of integrators
    in the open-loop transfer function using sympy.
    
    Args:
    num (list or sympy expression): Numerator coefficients or expression of the transfer function.
    den (list or sympy expression): Denominator coefficients or expression of the transfer function.
    
    Returns:
    dict: System type, description of steady-state error, and recommended steps.
    """
    # Define the Laplace variable 's'
    s = sp.symbols('s')
    
    # Create symbolic expressions for the transfer function's numerator and denominator
    num_expr = sum(coef * s**i for i, coef in enumerate(reversed(num)))
    den_expr = sum(coef * s**i for i, coef in enumerate(reversed(den)))
    
    # Solve for the poles (set the denominator equal to 0 and solve for s)
    poles = sp.solve(den_expr, s)
    
    # Count the number of integrators (poles at the origin, s = 0)
    num_integrators = sum(1 for pole in poles if abs(pole) < 1e-6)
    
    # Determine the system type based on the number of integrators
    system_type = num_integrators
    
    # Define the steady-state error expectations based on the system type
    if system_type == 0:
        step_error = "Finite error (proportional to 1/(1 + K))"
        ramp_error = "Infinite error"
        parabola_error = "Infinite error"
    elif system_type == 1:
        step_error = "Zero error"
        ramp_error = "Finite error (proportional to 1/K)"
        parabola_error = "Infinite error"
    elif system_type == 2:
        step_error = "Zero error"
        ramp_error = "Zero error"
        parabola_error = "Finite error"
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

def main1():
    # Example usage with a sample open-loop transfer function:
    # Define the numerator and denominator of the open-loop transfer function P(s)C(s) = K / (s(s+1)(s+2))
    K = 1
    num = [K]  # Numerator of the transfer function
    den = [1, 3, 2, 0]  # Denominator of the transfer function (with an integrator)

    # Calculate the system type and print the result
    system_info = calculate_system_type_by_array(num, den)
    print(system_info)

def calculate_system_type_by_denominator(den_expr: sp.core.mul.Mul):
    """
    Optimized system type calculation that directly checks for multiples of s in the denominator.
    """
    # Define the Laplace variable 's'
    s = sp.symbols('s')
    # print("tf:", den_expr)
    simplified_den = den_expr

    # Check if the denominator is a simple multiple of s (i.e., contains only powers of s)
    # simplified_den = sp.simplify(den_expr)
    # print("Simplified:", simplified_den)
    # Try to extract s**n if the denominator is purely a multiple of s
    # try:
    #     # Get the degree of s in the denominator, if it exists
    #     degree_of_s = sp.Poly(simplified_den, s).degree()
    # except:
    #     # If it's not a polynomial in s, fall back to solving
    #     degree_of_s = None

    # If we can extract the degree of s, we don't need to solve
    # if degree_of_s is not None:
    #     num_integrators = degree_of_s
    # else:
        # Solve for the poles (set the denominator equal to 0 and solve for s)
    # poles = sp.solve(simplified_den, s)
    # Step 1: Extract the denominator of the transfer function

    ###
    # denominator = simplified_den

    # # Step 2: Simplify the denominator (optional, if necessary)
    # simplified_den = sp.simplify(denominator)

    ###

    # Step 3: Use sp.roots() to get the poles and their multiplicities
    poles_with_multiplicities = sp.roots(simplified_den, s)

    # Step 4: Convert the dictionary of poles and multiplicities into a list that includes all multiples
    poles = []
    for pole, multiplicity in poles_with_multiplicities.items():
        poles.extend([pole] * multiplicity)  # Add the pole multiple times based on its multiplicity

    # Output the results
    # print(f"Simplified denominator: {simplified_den}")
    # print(f"Poles with multiplicities: {poles_with_multiplicities}")
    # print(f"All poles including multiplicities: {poles}")

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




import sympy as sp

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
        tracking_error_expr = 1 / ((s**q_cur) * tf)
        
        # Take the limit as s approaches 0 (final steady-state error)
        final_tracking_error = sp.limit(tracking_error_expr, s, 0)
        
        # Store the result in the dictionary with q_cur as the key
        tracking_errors[q_cur] = final_tracking_error
    
    return tracking_errors


def steady_state_error(tf: sp.core.mul.Mul, d: sp.core.mul.Mul):
    s = sp.symbols('s')
    return sp.limit(tf*d,s,0)





def main2():
    # Example usage with symbolic transfer function:
    s = sp.symbols('s')
    
    # Define the open-loop transfer function symbolically
    P = 1 / (s * (s + 1) * (s + 2))  # Plant transfer function
    C = 1  # Controller transfer function (simple proportional controller)

    transfer_function = P * C

    # Get the denominator of the transfer function
    den_expr = sp.denom(transfer_function)

    # Calculate the system type based on the denominator
    system_info = calculate_system_type_by_denominator(den_expr)
    print(system_info)

if __name__ == '__main__':
    # main1()  # Run the first main function for array-based input
    main2()  # Run the second main function for symbolic input (uncomment to use)



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
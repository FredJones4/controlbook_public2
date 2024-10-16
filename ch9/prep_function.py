import control as ctrl

def calculate_system_type(open_loop_tf):
    """
    Calculate the system type based on the number of integrators
    in the open-loop transfer function.
    
    Args:
    open_loop_tf (TransferFunction): Open-loop transfer function in s-domain.
    
    Returns:
    dict: System type, description of steady-state error, and recommended steps.
    """
    # Get the poles of the open-loop transfer function
    # Review of Control library Poles function: 
    poles = ctrl.poles(open_loop_tf)
    
    # Count the number of integrators (poles at the origin, s=0)
    # For those who need a review in Python List Comprenehsion:
    # https://www.w3schools.com/python/python_lists_comprehension.asp
    num_integrators = sum(1 for pole in poles if abs(pole) < 1e-6)
    
    # Determine the system type based on the number of integrators
    system_type = num_integrators
    
    # Define the steady-state error expectations based on the system type
    if system_type == 0:
        step_error = f"Finite error (proportional to 1/(1 + K))"
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

    # Output
    result = {
        "System Type": system_type,
        "Step Input Error": step_error,
        "Ramp Input Error": ramp_error,
        "Parabola Input Error": parabola_error,
        "Recommendation": recommendations.get(system_type, "No further improvements needed.")
    }
    
    return result
def main():
    # Example usage with a sample open-loop transfer function:
    # Define the open-loop transfer function P(s)C(s) = K / (s(s+1)(s+2))
    K = 1
    open_loop_tf = ctrl.TransferFunction([K], [1, 3, 2, 0])  # An example with an integrator

    # Calculate the system type and print the result
    system_info = calculate_system_type(open_loop_tf)
    print(system_info)

if __name__ == '__main__':
    main()

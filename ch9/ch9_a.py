from sympy_system_type import *
# Transfer function work
s = sp.symbols('s')
# other related symbols
k_d, k_p, k_i, ell, m, b = sp.symbols('k_d k_p k_i ell m b')
P = (3/m*ell**2) / (s*(s+3*b/(m*ell**2)))
C_pid = (k_d*s**2 + k_p*s + k_i)/s
C_pd = k_d*s + k_p
# Question A
# 
print(""""\033[1mQuestion A-1: When the controller for the single link robot arm is PD control, what is the
system type? Characterize the steady state error when the reference input
is a step, a ramp, and a parabola."\033[0m""")
# PD Form: System Type
transfer_function_pd = P*C_pd
den_expr = sp.denom(transfer_function_pd)
system_info_pd = calculate_system_type_by_denominator(den_expr)
print(system_info_pd)
# PD Form: Tracking Error
q_pd = system_info_pd['System Type']
tracking_error_pd = tracking_error(transfer_function_pd, q_pd)
print("Tracking Error for PD form: ", tracking_error_pd)
print(""""\033[1mQuestion A-2:  How does this change if you add an integrator?"\033[0m""")
# PID Form: System Type
transfer_function_pid = P*C_pid
den_expr = sp.denom(transfer_function_pid)
system_info_pid = calculate_system_type_by_denominator(den_expr)
print(system_info_pid)
print(system_info_pd)
# PID Form: Tracking Error
q_pid = system_info_pid['System Type']
tracking_error_pid = tracking_error(transfer_function_pid, q_pid)
print("Tracking error for PID form:", tracking_error_pid)
# Question B
# 
print(""""\033[1mQuestion B:  Consider the case where a constant disturbance acts at the input to the plant
(for example gravity in this case). What is the steady state error to a constant
input disturbance when the integrator is not present **and when it is present?"\033[0m""")
# PD Form
tf_d = s*P/(1+P*C_pd)
den_expr = sp.denom(tf_d)
q_d_pd = calculate_system_type_by_denominator(den_expr)['System Type']
print("System type:", q_d_pd)
A = sp.symbols('A')
D = A/(s**(q_d_pd+1))
sse_pd = steady_state_error(tf_d,D)
print("PD Form, Steady State error: ", sse_pd)
# PID Form




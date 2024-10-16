import math

# Part (a)
print("Part (a):")

# Given values
zeta = 0.7
tr = 2  # seconds

# Calculate natural frequency
wn = 2.2 / tr
print(f"Natural frequency (ωn) = {wn:.2f} rad/s")

# Calculate desired closed loop characteristic polynomial
s2_coeff = 1
s1_coeff = 2 * zeta * wn
s0_coeff = wn ** 2

print(f"Desired characteristic polynomial: s^2 + {s1_coeff:.2f}s + {s0_coeff:.2f}")

# Calculate kP and kD
kP = (s0_coeff - 0.6) / 0.2
kD = (s1_coeff - 0.1) / 0.2

print(f"Proportional gain (kP) = {kP:.2f}")
print(f"Derivative gain (kD) = {kD:.2f}")

# Part (b)
print("\nPart (b):")

# Given values
Fmax = 6  # N

# Calculate kP based on saturation constraint
kP_b = Fmax

print(f"New proportional gain (kP) = {kP_b:.2f}")

# Calculate kD to maintain damping ratio of 0.7
wn_b = math.sqrt(0.6 + 0.2 * kP_b)
kD_b = (2 * zeta * wn_b - 0.1) / 0.2

print(f"New derivative gain (kD) = {kD_b:.2f}")
print(f"New natural frequency (ωn) = {wn_b:.2f} rad/s")
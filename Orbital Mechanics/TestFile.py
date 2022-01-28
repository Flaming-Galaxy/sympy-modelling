import numpy as np
from TwoBody import two_body

# Set up initial conditions
mass = np.array([1e26, 1e26])

R1_0 = np.array([0, 0, 0])
R2_0 = np.array([3000, 0, 0])
V1_0 = np.array([10, 20, 30])
V2_0 = np.array([0, 40, 0])

y0 = np.concatenate((R1_0, R2_0, V1_0, V2_0))

tspan = np.array([0, 480])

answer = two_body(mass, y0, tspan)
print(answer)
import numpy as np
import numpy.linalg as linalg

# Initial conditions
R0 = np.array([8182.4, -6865.9, 0])
V0 = np.array([0.47572, 8.8116, 0])

# Gravitational parameter
G = 6.67259e-20 # Gravitational constant
mass = np.array([5.974e24, 1000])
mu = G*np.sum(mass)

# Other parameters
dtheta = 2*np.pi/3 # Radians
h = linalg.norm(np.cross(R0, V0)) # Angular momentum
Vr0 = np.dot(V0, R0)/linalg.norm(R0) # Initial radial velocity
R0 = linalg.norm(R0) # Initial radius
S = np.sin(dtheta)
C = np.cos(dtheta)

# Something
R = h**2/(mu*(1+(h**2/(mu*R0) - 1)*C - h*Vr0*S/mu))

# Compute Lagrange coefficients
f = 1 - mu*R*(1-C)/h**2
g = R*R0*S/h
fdot = (mu/h)*((1-C)/S)*((mu/h**2)*(1-C)-1/R0-1/R)
gdot = 1 - mu*R0*(1-C)/h**2

# Calculate position and velocity vector after change in true anomaly
R = f*R0 + g*V0
V = fdot*R0 + gdot*V0

# Return results
print("After {} radians change in true anomaly, the position vector is {} and the velocity vector is {}".format(dtheta, R, V))
import numpy as np
from scipy.integrate import solve_ivp
import numpy.linalg as linalg
import matplotlib.pyplot as plt

R0 = np.array([450+6378, 0, 0])
V0 = np.array([0, 7.650774, 0])

tspan = (0, 50*3600)

y0 = np.concatenate((R0, V0))

G = 6.67259e-20 # Gravitational constant
mass = np.array([5.974e24, 1000])
mu = G*(mass.sum())

def orbit(y0, tspan):
    
    def ode_func(t, y):
        # Get position vector of satellite
        Rs = y[:3]
        
        # Get velocity vector of satellite
        Vs = y[3:6]
        
        # Calculate acceleration vector of satellite
        r = linalg.norm(Rs)
        rho = 1.585e-12*np.exp(-(linalg.norm(Rs)-450)/62.2)
        C_D = 2.2
        A = 22400
        a_p = -rho*C_D*A*linalg.norm(Vs)/(2*mass[1])
        
        As = -mu*Rs/(r**3) + [a_p, 0, 0]
        
        # Return input for RKF4(5) method
        dydt = np.concatenate((Vs, As))
        return dydt
    
    sol = solve_ivp(ode_func, tspan, y0, method="RK45", max_step=5)
    return sol

sol = orbit(y0, tspan)

# Plot Earth as a wireframe sphere
R = 6378
fig = plt.figure("Orbit of 1000kg satellite")
ax = fig.add_subplot(projection='3d')

# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# x = R*np.cos(u)*np.sin(v)
# y = R*np.sin(u)*np.sin(v)
# z = R*np.cos(v)
# ax.plot_wireframe(x, y, z)

# Find position vector
Rx = sol.y[0]
Ry = sol.y[1]
Rz = sol.y[2]

# Plot satellite orbit
ax.plot(Rx, Ry, Rz)
plt.show()
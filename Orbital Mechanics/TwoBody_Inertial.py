# Numerical solution of the two body problem relative to an inertial frame

from scipy.integrate import solve_ivp
import numpy.linalg as linalg

def two_body(mass, y0, tspan):
    # Gravitational constant
    G = 6.67259e-20
    
    # Position vectors of m1 and m2
    R1 = y0[0:3]
    R2 = y0[3:6]
    
    # Velocity vectors of m1 and m2
    V1 = y0[6:9]
    V2 = y0[9:12]
    
    # Acceleration vectors of m1 and m2
    r = linalg.norm(R2-R1)
    A1 = G*mass[1]*(R2-R1)/(r**3)
    A2 = G*mass[0]*(R2-R1)/(r**3)
    
    # Integrate equations of motion
    dydt = [V1, V2, A1, A2]
    answer = solve_ivp(dydt, tspan, y0, method="RK45")
    
    return answer
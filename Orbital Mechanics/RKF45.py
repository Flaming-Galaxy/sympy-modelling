#Solve systems of first order differential equations
#Variable step size using RKF4(5) method

import numpy as np

def RKF45(ode_function, tspan, y0, tol=1e-8):
    # * Checking inputs
    # Checking ode_function is a vector
    ode_size = np.shape(ode_function)
    if not(ode_size[0] == 1 | ode_size[1] == 1):
        return "ode_function should be a vector"
    
    # Checking tspan is correct size (2,1)
    tspan_size = np.shape(tspan)
    if not(tspan_size == (2,1) | tspan_size == (1,2)):
        return "tspan should be a 2x1 vector of t0 and t_final"
    
    # Checking y0 is a vector
    y0_size = np.shape(y0)
    if not(y0_size[0] == 1 | y0_size[1] == 1):
        return "y0 should be a vector of initial conditions"
    
    # * Values for calculations
    a = np.array([0, 1/4, 3/8, 12/13, 1, 1/2])
    b = np.array([[0, 0, 0, 0, 0],
                  [1/4, 0, 0, 0, 0],
                  [3/32, 9/32, 0, 0, 0],
                  [1932/2197, -7200/2197, 7296/2197, 0, 0],
                  [439/216, -8, 3680/513, -845/4104, 0],
                  [-8/27, 2, -3544/2565, 1859/4104, -11/40]])
    c4 = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
    c5 = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
    
    # * Set up initial values
    t0 = tspan[0]
    t_final = tspan[1]
    t = t0
    y = y0
    h = (t_final - t0)/100
    
    # * Solve equations
    while (t < t_final):
        h_min = 16*np.finfo(np.float64).eps
        ti = t 
        yi = y 
        
        # Evaluate time derivatives at 6 points in the interval
        for i in range(6):
            t_inner = ti+a[i]*h
            y_inner = yi
            
            for j in range(5):
                y_inner = y_inner + h*b[i, j]*f[:,j]               
        
#Modelling a double pendulum system fixed at a point

from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point

#Create symbols
l_scalar = symbols('l')
theta_1_scalar, theta_2_scalar, x_scalar = dynamicsymbols('theta_1, theta_2, x')
x_scalar_dot = dynamicsymbols('x', 1)

#Create reference frames
N = ReferenceFrame('N')
A = ReferenceFrame('A')
B = ReferenceFrame('B')
A.orient(N, 'Axis', (theta_1_scalar, N.z))
B.orient(A, 'Axis', (theta_2_scalar, A.z))

#Angular velocities
N_w_A_vector = A.ang_vel_in(N)
N_w_B_vector = B.ang_vel_in(N)
A_w_B_vector = B.ang_vel_in(A)

#Angular accelerations
N_alpha_A_vector = A.ang_acc_in(N)
N_alpha_B_vector = B.ang_acc_in(N)
A_alpha_B_vector = B.ang_acc_in(A)

#Create points
O = Point('O')
Q = Point('Q')
P = Point('P')

#Define positions of points relative to each other
Q.set_pos(O, l_scalar*A.y)
P.set_pos(Q, x_scalar*B.y)

#Define vectors between points
r_OQ_vector = Q.pos_from(O)
r_QP_vector = P.pos_from(Q)
r_OP_vector = r_OQ_vector + r_QP_vector

#Set and store known velocities
O.set_vel(N, 0)
O.set_vel(A, 0)
N_v_O_vector = O.vel(N)
N_a_O_vector = O.acc(N)

Q.set_vel(A, 0)
Q.set_vel(B, 0)

P.set_vel(B, x_scalar_dot*B.y)
B_v_P_vector = P.vel(B)
B_a_P_vector = P.acc(B)
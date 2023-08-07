import numpy as np
import sympy as sp
from scipy.optimize import fsolve, root
import time
from numpy import pi, degrees, sin, cos

L1, L2, L3, L4, theta1 = (1, 2.5, 2.5, 2, 5/4*pi)
sin1, cos1 = (sin(theta1), cos(theta1))


def func_angles(vars):
    theta2, theta3 = vars
    return [
        L1*sin1 + L2*sin(theta2) + L3*sin(theta3),
        L1*cos1 + L2*cos(theta2) + L3*cos(theta3) - L4
    ]

def jac_angles(vars):
    theta2, theta3 = vars
    return [
        [L2*cos(theta2), L3*cos(theta3)],
        [-L2*sin(theta2), -L3*sin(theta3)]
    ]

def func_sincos(vars):
    sin2, cos2, sin3, cos3 = vars
    return [
        L1*sin1 + L2*sin2 + L3*sin3,
        L1*cos1 + L2*cos2 + L3*cos3 - L4, 
        sin2**2 + cos2**2 - 1,
        sin3**2 + cos3**2 - 1,        
    ]

def jac_sincos(vars):
    sin2, cos2, sin3, cos3 = vars
    return [
        [L2, 0, L3, 0],
        [0, L2, 0, L3],
        [2*sin2, 2*cos2, 0, 0],
        [0, 0, 2*sin3, 2*cos3]
    ]


# Generate random initial guesses
guess = np.array([1, 5])
guess4 = np.array([sin(guess[0]), cos(guess[0]), sin(guess[1]), cos(guess[1])])



time1, time2 = (0, 0)

n = 1000
for i in range(n):


    start_time = time.time()
    solution1 = root(func_sincos, jac=jac_sincos, x0=guess4, method='hybr')
    # solution1 = np.mod(solution1.x, 2*pi)
    solution1 = solution1.x
    time1 += time.time() - start_time


    start_time = time.time()
    x0 = guess4
    tolerance = 1e-6
    max_iterations = 100
    for iteration in range(max_iterations):
        sin2, cos2, sin3, cos3 = x0
        h = [
            (L2*cos3*(cos2**2 + sin2**2 - 1) + L3*cos2*(cos3**2 + sin3**2 - 1) - 2*cos2*(cos3*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin3*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L2*(cos2*sin3 - cos3*sin2)),
            (-L2*sin3*(cos2**2 + sin2**2 - 1) - L3*sin2*(cos3**2 + sin3**2 - 1) + 2*sin2*(cos3*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin3*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L2*(cos2*sin3 - cos3*sin2)),
            (-L2*cos3*(cos2**2 + sin2**2 - 1) - L3*cos2*(cos3**2 + sin3**2 - 1) + 2*cos3*(cos2*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin2*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L3*(cos2*sin3 - cos3*sin2)),
            (L2*sin3*(cos2**2 + sin2**2 - 1) + L3*sin2*(cos3**2 + sin3**2 - 1) - 2*sin3*(cos2*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin2*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L3*(cos2*sin3 - cos3*sin2))
        ]
        x_new = x0 + h

        # print(x_new)

        if np.linalg.norm(x_new - x0) < tolerance:
            break

        x0 = x_new
    solution2 = x0
    time2 += time.time() - start_time





print("1:    ", solution1, time1*1000/n)
print("2:    ", solution2, time2*1000/n)



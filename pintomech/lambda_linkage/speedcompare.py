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
guess2 = np.array([1, 5])
guess4 = np.array([sin(guess2[0]), cos(guess2[0]), sin(guess2[1]), cos(guess2[1])])


time1, time2 = (0, 0)

def func_sincos(vars):
    sin2, cos2, sin3, cos3 = vars
    sin1 = sin(theta1)
    cos1 = cos(theta1)
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

ntests = 4
times = np.zeros(ntests)

n = 1000
for i in range(n):


    start_time = time.time()
    sol0 = root(func_sincos, jac=jac_sincos, x0=guess4, method='hybr').x
    times[0] += time.time() - start_time

    start_time = time.time()
    sol1 = root(func_angles, jac=jac_angles, x0=guess2).x
    times[1] += time.time() - start_time

    start_time = time.time()
    tolerance = 1e-6
    max_iterations = 100
    x0 = guess2
    for iteration in range(max_iterations):
        theta2, theta3 = x0
        h = [
            (L1*cos(theta1-theta3) + L2*cos(theta2-theta3) + L3 - L4*cos(theta3)) / (L2*(theta2-theta3)),
            (-L1*cos(theta1-theta2) - L2 - L3*cos(theta2-theta3) + L4*cos(theta2)) / (L3*(theta2-theta3))
        ]
        x_new = x0 + h
        if np.linalg.norm(x_new - x0) < tolerance:
            break
        x0 = x_new
    sol2 = x0
    times[2] += time.time() - start_time

    start_time = time.time()
    tolerance = 1e-6
    max_iterations = 100
    x0 = guess4
    for iteration in range(max_iterations):
        sin2, cos2, sin3, cos3 = x0
        h = [
            (L2*cos3*(cos2**2 + sin2**2 - 1) + L3*cos2*(cos3**2 + sin3**2 - 1) - 2*cos2*(cos3*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin3*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L2*(cos2*sin3 - cos3*sin2)),
            (-L2*sin3*(cos2**2 + sin2**2 - 1) - L3*sin2*(cos3**2 + sin3**2 - 1) + 2*sin2*(cos3*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin3*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L2*(cos2*sin3 - cos3*sin2)),
            (-L2*cos3*(cos2**2 + sin2**2 - 1) - L3*cos2*(cos3**2 + sin3**2 - 1) + 2*cos3*(cos2*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin2*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L3*(cos2*sin3 - cos3*sin2)),
            (L2*sin3*(cos2**2 + sin2**2 - 1) + L3*sin2*(cos3**2 + sin3**2 - 1) - 2*sin3*(cos2*(L1*cos1 + L2*cos2 + L3*cos3 - L4) + sin2*(L1*sin1 + L2*sin2 + L3*sin3)))/(2*L3*(cos2*sin3 - cos3*sin2))
        ]
        x_new = x0 + h
        if np.linalg.norm(x_new - x0) < tolerance:
            break
        x0 = x_new
    sol3 = x0
    times[3] += time.time() - start_time



for i in range(ntests):
    print(i, round(times[i], 5))
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi, arccos, sqrt, arctan

# Define the coordinates of the end effector
x = 0.5  # x-coordinate
y = -0.8  # y-coordinate

# Define the link lengths of the robot arm
a1 = 1
a2 = 1

# Calculate theta2 using inverse kinematics
D = (x**2 + y**2 - a1**2 - a2**2) / (2 * a1 * a2)
theta2 = np.arctan2(-np.sqrt(1 - D**2), D)
theta2_jack = pi - arccos((a1**2 + a2**2 - x**2 - y**2) / (2*a1*a2))
theta1_jack = 2*pi - arccos((a1**2 - a2**2 + x**2 + y**2) / (2*a1*sqrt(x**2 + y**2))) + arctan(y/x)

# Calculate theta1 using inverse kinematics
theta1 = np.arctan2(y, x) - np.arctan2(a2 * np.sin(theta2), a1 + a2 * np.cos(theta2))

# Plot the robot arm
plt.figure()
plt.plot([0, a1*np.cos(theta1), a1*np.cos(theta1) + a2*np.cos(theta1 + theta2)], [0, a1*np.sin(theta1), a1*np.sin(theta1) + a2*np.sin(theta1 + theta2)], 'bo-')


plt.plot([0, a1*np.cos(theta1_jack), a1*np.cos(theta1_jack) + a2*np.cos(theta1_jack + theta2_jack)], [0, a1*np.sin(theta1_jack), a1*np.sin(theta1_jack) + a2*np.sin(theta1_jack + theta2_jack)], 'go-')
plt.plot(0, 0, 'o', color='orange')  # Mark the origin as orange
plt.text(0, 0, f'Theta1: {np.degrees(theta1):.2f} deg', fontsize=6, color='red')  # Display theta1
plt.text(1.5, 2, f'Theta2: {np.degrees(theta2):.2f} deg', fontsize=6, color='blue')  # Display theta2

plt.text(-0.5, 0.5, f'Theta1_jack: {np.degrees(theta1_jack):.2f} deg', fontsize=6, color='red')  # Display theta1
plt.text(-1.5, 2, f'Theta2_jack: {np.degrees(theta2_jack):.2f} deg', fontsize=6, color='blue')  # Display theta2
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
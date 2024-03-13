import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, arccos, sqrt, arctan


# Define the link lengths of the robot arm
a1 = 1
a2 = 1

def ik(x, y): #input is coords of end effector
    D = (x**2 + y**2 - a1**2 - a2**2) / (2 * a1 * a2)
    theta2 = np.arctan2(-np.sqrt(1 - D**2), D)
    theta1 = np.arctan2(y, x) - np.arctan2(a2 * np.sin(theta2), a1 + a2 * np.cos(theta2))
    return [theta1, theta2]

def ik_jack(x, y): #input is coords of end effector
    theta1 = 2*pi - arccos((a1**2 - a2**2 + x**2 + y**2) / (2*a1*sqrt(x**2 + y**2))) + np.arctan2(y,x)
    theta2 = pi - arccos((a1**2 + a2**2 - x**2 - y**2) / (2*a1*a2))
    return [theta1, theta2]


class DraggableCircle:
    def __init__(self, ax):
        self.ax = ax

        self.x = 0.5 #starting end effector pos
        self.y = -0.8

        self.circle = plt.Circle((self.x, self.y), 0.1, fc='b', alpha=0.5)  #center, radius, facecolor, alpha
        self.ax.add_patch(self.circle)

        [theta1, theta2] = ik(self.x, self.y)
        self.arm1, = ax.plot(
            [0, a1*np.cos(theta1), a1*np.cos(theta1) + a2*np.cos(theta1 + theta2)], #x coords of arm
            [0, a1*np.sin(theta1), a1*np.sin(theta1) + a2*np.sin(theta1 + theta2)], #y coords of arm
            'go-'
        )
        [theta1, theta2] = ik_jack(self.x, self.y)
        self.arm1_jack, = ax.plot(
            [0, a1*np.cos(theta1), a1*np.cos(theta1) + a2*np.cos(theta1 + theta2)], #x coords of arm
            [0, a1*np.sin(theta1), a1*np.sin(theta1) + a2*np.sin(theta1 + theta2)], #y coords of arm
            'bo-'
        )

        self.cid_press = self.circle.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.circle.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.circle.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.dragging = False

    def on_press(self, event):
        if event.inaxes != self.circle.axes:
            return
        contains, _ = self.circle.contains(event)
        if not contains:
            return
        self.dragging = True
        self.offset = (self.circle.center[0] - event.xdata, self.circle.center[1] - event.ydata)

    def on_release(self, event):
        self.dragging = False

    def on_motion(self, event):
        if not self.dragging or event.xdata is None:
            return
        
        (self.x, self.y) = (event.xdata + self.offset[0], event.ydata + self.offset[1])
        self.circle.center = (self.x, self.y)

        try:
            [theta1, theta2] = ik(self.x, self.y)
            self.arm1.set_xdata([0, a1*np.cos(theta1), a1*np.cos(theta1) + a2*np.cos(theta1 + theta2)]) #x coords of arm
            self.arm1.set_ydata([0, a1*np.sin(theta1), a1*np.sin(theta1) + a2*np.sin(theta1 + theta2)]) #y coords of arm)

            [theta1, theta2] = ik_jack(self.x, self.y)
            self.arm1_jack.set_xdata([0, a1*np.cos(theta1), a1*np.cos(theta1) + a2*np.cos(theta1 + theta2)]) #x coords of arm
            self.arm1_jack.set_ydata([0, a1*np.sin(theta1), a1*np.sin(theta1) + a2*np.sin(theta1 + theta2)]) #y coords of arm)
        except:
            print('out of range')
            return

        self.circle.figure.canvas.draw()


plt.figure()
plt.title("Drag the purple")
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
dc = DraggableCircle(ax)
plt.show()
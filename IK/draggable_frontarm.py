import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, arccos, sqrt, arctan2, sin, cos
import warnings


# link lengths of front arm (mm), refer to Jack's IK sketch
a1 = 23.785
a2 = 72.253
a3 = 35.190
a4 = 47.276
a5 = 29.621
a6 = 25.939
a7 = 20.200

def HTM(theta, x, y):
    return np.array([
        [cos(theta), -sin(theta), x],
        [sin(theta), cos(theta), y],
        [0, 0, 1]
    ])
    

def ik(x, y): #input is coords of end effector
    theta2 = pi - arccos((a1**2 + a2**2 - x**2 - y**2) / (2*a1*a2))
    theta1 = 2*pi - arccos((a1**2 - a2**2 + x**2 + y**2) / (2*a1*sqrt(x**2 + y**2))) + arctan2(y,x)
    # beta = arccos((a2**2 + a3**2 - a4**2) / (2*a2*a3))
    # O2 = np.array([a1*cos(theta1), a1*sin(theta1)])

    # H_2_to_4 = np.array([
    #     [cos(theta1+beta), -sin(theta1+beta), a3*cos(theta1+beta)],
    #     [sin(theta1+beta), cos(theta1+beta), a3*sin(theta1+beta)],
    #     [0, 0, 1]
    # ])
    # O4 = O2 + H_2_to_4@O2
    

    # theta5 = pi + arccos((a5**2 + a6**2 - x4**2 - y4**2) / (-2*a5*a6))
    # theta6 = arccos((a6**2 + x4**2 + y4**2 - a5**2) / (2*a6*sqrt(x4**2 + y4**2))) + arctan2(y4,x4)

    return [theta1, theta2]
    return


class DraggableCircle:
    def __init__(self, ax):
        self.ax = ax

        self.x = 56.073 #starting end effector pos
        self.y = 19.815

        self.circle = plt.Circle((self.x, self.y), 10, fc='b', alpha=0.5)  #center, radius, facecolor, alpha
        self.ax.add_patch(self.circle)

        [theta1, theta2] = ik(self.x, self.y)
        self.arm1, = ax.plot(
            [0]*6, #x coords of arm
            [0]*6, #y coords of arm
            'go-',
            markersize=3
        )

        self.O4 = ax.scatter(150,150, c='orange')

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


        [x,y] = [self.x, self.y]
        # print(x, y)

        try:

            H0 = HTM(0, x=0, y=0)

            Heff = np.linalg.inv(H0) @ np.array([x,y,1])
            [x,y] = [Heff[0], Heff[1]]

            theta2 = pi - arccos((a1**2 + a2**2 - x**2 - y**2) / (2*a1*a2))
            theta1 = (2*pi - arccos((a1**2 - a2**2 + x**2 + y**2) / (2*a1*sqrt(x**2 + y**2))) + arctan2(y,x)) % (2*pi)
            
            beta = arccos((a2**2 + a3**2 - a4**2) / (2*a2*a3))

            H1 = H0 @ HTM(theta1, x=0, y=0)
            H2 = H1 @ HTM(theta2, x=a1, y=0)
            H3 = H2 @ HTM(0, x=a2, y=0)
            H4 = H2 @ HTM(0, x=a3*cos(beta), y=a3*sin(beta))


            H4_shifted = np.linalg.inv(H0 @ HTM(0, 0, a7)) @ H4 #get H4 relative to second servo
            [x4,y4] = [H4_shifted[0,2], H4_shifted[1,2]]

            D = (x4**2 + y4**2 - a6**2 - a5**2) / (2 * a6 * a5)
            theta5 = np.arctan2(-np.sqrt(1 - D**2), D)
            theta6 = np.arctan2(y4, x4) - np.arctan2(a5 * np.sin(theta5), a6 + a5 * np.cos(theta5))

            # theta5 = pi + arccos((a5**2 + a6**2 - x4**2 - y4**2) / -(2*a5*a6))
            # theta6 = arccos((a6**2 - a5**2 + x4**2 + y4**2) / (2*a6*sqrt(x4**2 + y4**2))) + arctan2(y4,x4)

            H6 = H0 @ HTM(theta6, x=0, y=a7)
            H5 = H6 @ HTM(theta5, x=a6, y=0)
            # H6 = H0 @ HTM(pi/6, x=0, y=a7)
            # H5 = H6 @ HTM(pi/2, x=a6, y=0)
            H4_check = H5 @ HTM(0, x=a5, y=0)

            self.arm1.set_xdata([H1[0,2], H2[0,2], H3[0,2], H4[0,2], H2[0,2], H4[0,2], H5[0,2], H6[0,2], H5[0,2], H4_check[0,2]]) #x coords of arm
            self.arm1.set_ydata([H1[1,2], H2[1,2], H3[1,2], H4[1,2], H2[1,2], H4[1,2], H5[1,2], H6[1,2], H5[1,2], H4_check[1,2]]) #y coords of arm

            # self.O4.set_offsets([x4, y4])
            self.ax.set_title(f'{np.round([theta1, theta6], 3)}')

            
        except Warning as w:
            print('unsolveable')

        self.circle.figure.canvas.draw()

warnings.catch_warnings()
warnings.simplefilter("error")

plt.figure()
plt.title("Drag the purple")
ax = plt.gca()

outer_circle = plt.Circle((0, 0), a1+a2, color='blue', alpha=0.5)
inner_circle = plt.Circle((0, 0), a2-a1, color='white', alpha=1)
ax.add_artist(outer_circle)
ax.add_artist(inner_circle)

ax.set_aspect('equal', adjustable='box')
ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
dc = DraggableCircle(ax)

#workspace




plt.show()
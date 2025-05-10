import matplotlib.pyplot as plt

class DraggableCircle:
    def __init__(self, ax):
        self.ax = ax
        self.circle = plt.Circle((0, 0), 0.1, fc='b', alpha=0.5)  # Initial circle at origin
        self.ax.add_patch(self.circle)
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
        if not self.dragging:
            return
        self.circle.center = (event.xdata + self.offset[0], event.ydata + self.offset[1])
        self.circle.figure.canvas.draw()

plt.figure()
ax = plt.gca()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
dc = DraggableCircle(ax)
plt.show()
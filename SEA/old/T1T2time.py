import numpy as np
import matplotlib.pyplot as plt

def RK4(f, x, u, dt):
    # Runge-Kutta 4 integration
    k1,log = f(x, u)
    k2,_ = f(x + (dt/2)*k1, u)
    k3,_ = f(x + (dt/2)*k2, u)
    k4,_ = f(x + dt*k3, u)
    return x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4), log

t_max = 0.2 #200ms
dt = 0.001


t = 0;
t_seg = [0, 0.05, 0.10, 0.15, 0.2]
seg = 0;
seg_width = t_seg[1] - t_seg[0]

T1_prof = [1, 1.5, 3.2, 4.5, 5]


plt.figure()

for t in np.arange(0, t_max+dt, step=dt):
    if(t > t_seg[seg+1]):
        seg += 1
        seg_width = t_seg[seg+1] - t_seg[seg]

    interp = (t - t_seg[seg]) / seg_width

    T1 = interp*T1_prof[seg+1] + (1-interp)*T1_prof[seg]
    
    print(round(t, 3), T1)

    plt.plot(t, T1, 'ro', markersize=1)

plt.scatter(t_seg, T1_prof)

plt.show()


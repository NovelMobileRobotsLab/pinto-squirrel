import numpy as np
import matplotlib.pyplot as plt

def RK4(f, x, u, dt):
    # Runge-Kutta 4 integration
    k1,log = f(x, u)
    k2,_ = f(x + (dt/2)*k1, u)
    k3,_ = f(x + (dt/2)*k2, u)
    k4,_ = f(x + dt*k3, u)
    return x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4), log

def T1(x_a):
    return x_a
def dT1dxa(x_a):
    return 1
def T2_inv(x_l):
    return x_l + 0.1*x_l**2
def dT2invdxl(x_l):
    return 1 + 0.2*x_l


def sea_model(x, u):
    u_a, u_l = u
    x_a, x_a_dot, x_l, x_l_dot = x

    k = 0.69 #spring constant for xlf=15
    # k = 2
    # k = 0.28 #spring constant for xlf=25
    m_a = 0.1 #mass of the actuator before spring
    m_l = 0.5 #mass of the load after spring 
    tau_stall = 10
    v_a_max = 20

    tau_a = -tau_stall/v_a_max*x_a_dot + u_a*tau_stall
    x_a_ddot = (tau_a + k*(T2_inv(x_l) - T1(x_a))*dT1dxa(x_a)) / m_a
    x_l_ddot = (0 - k*(T2_inv(x_l) - T1(x_a))*dT2invdxl(x_l)) / m_l
    # x_a_ddot = (tau_a + k*(x_l-x_a)) / m_a
    # x_l_ddot = (tau_l - k*(x_l-x_a)) / m_l

    x_dot = np.array([x_a_dot, x_a_ddot, x_l_dot, x_l_ddot])
    log = [tau_a]
    return x_dot, log



disp = True
# disp = False

if(disp):
    import pygame
    pygame.init()
    screen = pygame.display.set_mode((1280, 720))
    clock = pygame.time.Clock()
    width = screen.get_width()
    height = screen.get_height()
    xrange = [-5, 40]
    pixels_per_meter = width / (xrange[1] - xrange[0])

    def to_disp_pos(x, y=0):
        return ((x-(xrange[1]+xrange[0])/2)*pixels_per_meter+width/2, -y*screen.get_height() + height/2)

running = True
dt = 1/60


t = 0
x = [0, 0, 0, 0]
u = [0, 0]

xs = []
logs = []
times = []

while running:
    
    if(x[2] > 15):
        running = False

    if(t > 0.5):
        u = [1, 0]

    x, log = RK4(sea_model, x, u, dt)

    t += dt
    times.append(t)
    xs.append(x)
    logs.append(log)

    if(disp):
        #check for exit
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        screen.fill("white")
        spring_offset = 0
        pygame.draw.line(screen, 'gray', (0,0), (100,100))
        pygame.draw.circle(screen, "red", to_disp_pos(x[0], 0), 1*pixels_per_meter)
        pygame.draw.circle(screen, "blue", to_disp_pos(x[2]+spring_offset, 0), 1*pixels_per_meter)
        pygame.draw.line(screen, "orange", to_disp_pos(x[0], 0), to_disp_pos(x[2]+spring_offset, 0), width=round(50/(abs(x[0]-x[2]-spring_offset)+1)))
        pygame.display.flip()
        dt = clock.tick(60) / 1000
        
if(disp):
    pygame.quit()

xs = np.array(xs)
logs = np.array(logs)

print(f"time taken: {t}")
plt.plot(times, xs[:,0], label='x_a')
plt.plot(times, xs[:,1], label='motor speed x_a_dot')
plt.plot(times, xs[:,2], label='x_l')
plt.plot(times, xs[:,3], label='x_l_dot')
plt.plot(times, logs[:,0], label='motor torque')
plt.plot(times, logs[:,0]*xs[:,1], label='power τ*ω')
plt.title(f"end load vel: {np.round(xs[-1,3], 2)}")
plt.legend()
plt.show()
from casadi import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc
import os
import datetime
from tqdm import tqdm
import itertools

topN = 5

tau_max=0.05
omega_max=3000
nsteps = 20
k_test = 2500
n_successes = 10 #takes about 4 seconds
n_max_trials = 20


m_l_range = np.linspace(0.100, 1.000, 5)
t_total_range = np.linspace(0.010, 0.300, 5)
leg_max_range = np.linspace(0.050, 0.2, 5)


permutations = list(itertools.product(m_l_range, t_total_range, leg_max_range))
permutations = np.array(permutations)

filedir = os.path.dirname(os.path.abspath(__file__))
timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
exp_path = f"{filedir}/exp_{timestamp}"
os.makedirs(exp_path, exist_ok=True)

for i in tqdm(range(len(permutations))):
    m_l = permutations[i,0]
    t_total = permutations[i,1]
    leg_max = permutations[i,2]
    # print(m_l, t_total, leg_max)
    param_str = f"k{round(k_test)}_m{round(m_l*1000)}_t{round(t_total*1000)}_L{round(leg_max*1000)}"

    t = t_total/nsteps #should be named dt

    ''' Function for one analytical step '''
    def get_Fstep():
        k = SX.sym('k')
        X0 = SX.sym('X0', 5) #states: x_a, x_l, v_l, T1_last, T2_last
        x_a0 = X0[0]
        x_l0 = X0[1]
        v_l0 = X0[2]
        T1_0 = X0[3]
        T2_0 = X0[4]
        T = SX.sym('T', 2) #transmission ratio for this step
        T1 = T[0]
        T2 = T[1]
        x_l0 = T2*(x_a0*(T1-T1_0) + x_l0/T2_0) #solve for x_l to make spring compression continuous

        x0 = 1/T2
        x1 = T2**2
        x2 = omega_max**2
        x3 = x1*x2
        x4 = k**2
        x5 = m_l**2
        x6 = T1**4*x5
        x7 = x4*x6
        x9 = 1/tau_max
        x10 = (1/2)*x9
        x12 = T1**2
        x13 = omega_max*t
        x14 = k*x10*x12*x13
        x15 = exp(x14)
        x16 = tau_max*x15
        x17 = x16*x_a0
        x18 = T1*omega_max
        x19 = T2*m_l*v_l0*x18
        x20 = m_l*x12
        x21 = x20*x3
        x22 = k*m_l
        x23 = 2*tau_max
        x24 = x23*x_l0
        x25 = x12*x22
        x26 = x23*x_a0
        x27 = T2*omega_max
        x28 = T2**3
        x29 = omega_max**3
        x30 = T1**3
        x32 = exp(-x14)
        x33 = x32*x9
        x34 = x1*x18
        x35 = T2*v_l0 - x34
        x36 = 2*k
        x37 = T1*T2
        x38 = omega_max*x1
        x39 = 2*v_l0*x25*x38
        x40 = x2*x28
        x41 = x22*x30*x40
        x42 = 2*x41
        x43 = tau_max*x36
        x44 = -x37*x43*x_a0 + x43*x_l0
        r = 4*k*m_l*tau_max**2 - x3*x7
        R_ad_re = if_else(r>0, sqrt(r), 0)
        R_ad_im = if_else(r>0, 0, sqrt(-r))
        x11_re = t*x0*x10*R_ad_re/m_l
        x11_im = t*x0*x10*R_ad_im/m_l
        C_os = cos(x11_re)*cosh(x11_im)
        S_in_re = sin(x11_re)*cosh(x11_im)
        S_in_im = cos(x11_re)*sinh(x11_im)
        SdR = (S_in_re*R_ad_re + S_in_im*R_ad_im) / (R_ad_re**2 + R_ad_im**2) # S_in/R_ad should be purely real (aka x31)
        SmR = S_in_re*R_ad_re - S_in_im*R_ad_im  # S_in*R_ad should be purely real
        x_a = x33*(C_os*(-x19 + x21) + x13*x16 + x15*x19 - x15*x21 + x17 + SdR*(-k*v_l0*x3*x30*x5 + k*x28*x29*x6 + x18*x22*x24 - x25*x26*x27))
        x_l = (1/2)*x33*(C_os*(-x39 + x42 + x44) + SmR*x35 + T2*t*x16*x18*x36 + x15*x39 - x15*x42 + x17*x36*x37 + SdR*(T1**5*T2**4*x29*x4*x5 - m_l*x26*x30*x38*x4 - v_l0*x40*x7 + x20*x24*x27*x4))/k
        v_l = x0*x32*(C_os*x35 + x15*x34 + SdR*(k*m_l*omega_max*v_l0*x1*x12 - x41 - x44))
        X_new = SX.sym('X_new', 5)
        X_new[0] = x_a
        X_new[1] = x_l
        X_new[2] = v_l
        X_new[3] = T1
        X_new[4] = T2

        Fstep = Function("step", [X0, T, k], [X_new])
        return Fstep
    Fstep_n = get_Fstep().mapaccum(nsteps)

    def opti_steps(nsteps, t, X_init, T_guess, k_guess, T1_bound, T2_bound, k_bound, leg_max):
        ''' Optimize T1 and T2 over multiple steps '''
        opti = Opti()

        Ts = opti.variable(2, nsteps) #column vectors of [T1, T2], concat horizontally
        k_var = opti.variable() #column vectors of [T1, T2], concat horizontally

        XT_init = vertcat(X_init, Ts[:,0])
        X_n = Fstep_n(XT_init, Ts, k_var)
        # X_n = Fstep_n(np.concatenate((X_init, T_guess[:,0])), Ts, k_var)

        opti.minimize(-X_n[2,-1])

        opti.set_initial(Ts[0,:], T_guess[0,0])
        opti.set_initial(Ts[1,:], T_guess[1,0])
        opti.set_initial(k_var, k_guess)

        opti.subject_to(opti.bounded(T1_bound[0], Ts[0,:], T1_bound[1]))
        opti.subject_to(opti.bounded(T2_bound[0], Ts[1,:], T2_bound[1]))
        # for i in range(nsteps-1):
        #     opti.subject_to(opti.bounded(-0.0001, Ts[0,i+1] - Ts[0,i], 0.0001))
        #     opti.subject_to(opti.bounded(-3, Ts[1,i+1] - Ts[1,i], 3))
        opti.subject_to(sum2(X_n[2,:])*t < leg_max) #limit the maximum leg extension using the integral of leg vel (sum across t)
        opti.subject_to(X_n[2,:] >= 0) #leg vel cannot be negative
        opti.subject_to(opti.bounded(k_bound[0], k_var, k_bound[1]))

        # opti.callback(lambda i: log.append(opti.debug.value(X_n))) #for logging each iteration

        opts = {'ipopt.print_level':0, 'print_time':0}
        opti.solver('ipopt', opts)
        sol = opti.solve()
        X_opt = np.array(sol.value(X_n))

        X_opt_first_col = np.append(X_init, [X_opt[3,0], X_opt[4,0]]).reshape(-1,1)
        X_opt = np.hstack([X_opt_first_col, X_opt]) #insert data at t=0
        k_opt = sol.value(k_var)

        return sol, X_opt, k_opt


    ''' Do trials '''

    lhc = qmc.LatinHypercube(d=2*nsteps+1)
    lhc_samp = lhc.random(n=n_max_trials)
    T1_bound = [1e-6, 0.002]
    T2_bound = [1e-6, 100]
    k_bound = [k_test-0.4, k_test+0.4]
    X_init = [0,0,0]

    param_path = f'{exp_path}/{param_str}'
    os.makedirs(param_path, exist_ok=True)

    successes = 0
    fails = 0
    i = 0
    while successes < n_successes and i < n_max_trials:
        # print(f"[{i}]", end='')
        try:
            T1_guess = lhc_samp[i,:nsteps]*(T1_bound[1] - T1_bound[0]) + T1_bound[0]
            T2_guess = lhc_samp[i,nsteps:2*nsteps]*(T2_bound[1] - T2_bound[0]) + T2_bound[0]
            k_guess = lhc_samp[i,2*nsteps]*(k_bound[1] - k_bound[0]) + k_bound[0]
            T_guess = np.array([T1_guess, T2_guess])

            sol, X_opt, k_opt = opti_steps(nsteps, t, X_init, T_guess, k_guess, T1_bound, T2_bound, k_bound, leg_max)
            end_v = X_opt[2,-1]

            np.savetxt(f'{param_path}/Xopt_v{round(end_v*1000)}_{param_str}.txt', X_opt)
            successes += 1
        except Exception as e:
            fails += 1
            # print(e)
        i += 1

    # print(f"successes: {successes}, fails: {fails}")

    #plot trials

    axs = [0,0]
    fig, axs[0] = plt.subplots(1,1, figsize=(10,6))
    ax0_twin = axs[0].twinx()

    v_max = 0
    for file in sorted(os.listdir(param_path))[-topN:]:
        if(not file.endswith(".txt")):
            continue

        X_optload = np.loadtxt(f'{param_path}/{file}')
        file = file[:-4] #remove '.txt'
        v_file = round(int(file.split("_")[1].split("v")[1]) * 0.001, 3)
        
        if(v_file > v_max):
            v_max = v_file
            X_best = X_optload

        ts = np.linspace(0, t_total, nsteps+1)
        axs[0].step(ts, X_optload[3], 'r', alpha=0.2)
        ax0_twin.step(ts, X_optload[4], 'b', alpha=0.2)

    axs[0].step(ts, X_best[3], 'r--', label=f"best", alpha=0.8, linewidth=2)
    ax0_twin.step(ts, X_best[4], 'b--', label=f"best", alpha=0.8, linewidth=2)
    ax0_twin.set_ylim(0, 20)
    axs[0].set_xlabel('time (seconds)')
    axs[0].set_ylabel('T1', color='red')
    ax0_twin.set_ylabel('T2', color='blue')

    plt.title(f"Top {topN} transmission profiles, best in dashed\nk={k_test}, m={m_l}, t={t_total}, L={leg_max}, v_best={v_max}")

    os.makedirs(f"{exp_path}/figs", exist_ok=True)
    plt.savefig(f'{exp_path}/figs/{param_str}.png', dpi=300)
    # plt.show()
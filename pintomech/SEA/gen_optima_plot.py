from casadi import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc
import os
import datetime
from tqdm import tqdm
import itertools

topN = 5
nsteps = 20

tau_max=0.05
omega_max=3000

filedir = os.path.dirname(os.path.abspath(__file__))
exp_name = "exp_2023-10-29-12-10-02"
exp_path = f"{filedir}/{exp_name}"


# param_folders = os.listdir(exp_path)
param_folders = [os.path.join(f"{exp_path}/params", folder) for folder in sorted(os.listdir(f"{exp_path}/params"))]
for param_path in param_folders:

    param_str = os.path.basename(param_path)
    print(param_str)

    if(not param_str.startswith('k')):
        continue

    k_test = int(param_str.split("_")[0].split("k")[1])
    m_l = round(int(param_str.split("_")[1].split("m")[1])*0.001, 3)
    t_total = round(int(param_str.split("_")[2].split("t")[1])*0.001, 3)
    leg_max = round(int(param_str.split("_")[3].split("L")[1])*0.001, 3)

    axs = [0,0]
    fig, axs[0] = plt.subplots(1,1, figsize=(10,6))
    # ax0_twin = axs[0].twinx()

    v_max = 0

    for file in sorted(os.listdir(param_path))[-topN:]:
        if(not file.endswith(".txt")):
            continue

        X_optload = np.loadtxt(f'{param_path}/{file}')
        file = file[:-4] #remove '.txt'
        v_file = round(int(file.split("_")[1].split("v")[1]) * 0.001, 3)
        
        

        ts = np.linspace(0, t_total, nsteps+1)
        t = t_total/nsteps

        x_as = X_optload[0,:]
        x_ls_adj = X_optload[1,:]
        v_ls = X_optload[2,:]
        T1s = X_optload[3,:]
        T2s = X_optload[4,:]

        v_file = v_ls[-1]
        if(v_file > v_max):
            v_max = v_file
            X_best = X_optload

        KE = (1/2)*m_l*v_ls**2
        PE = (1/2)*k_test*(x_as*T1s - x_ls_adj/T2s)**2

        # axs[0].step(ts, X_optload[2], 'r', alpha=0.2)

        axs[0].plot(ts, KE, 'blue', alpha=0.3)
        # axs[0].plot(ts, PE, 'orange', alpha=0.3)
        # axs[0].plot(ts, E_motor, 'red', alpha=0.2)
        axs[0].plot(ts, PE, 'orange', alpha=0.3)

    x_as = X_best[0,:]
    x_ls_adj = X_best[1,:]
    v_ls = X_best[2,:]
    T1s = X_best[3,:]
    T2s = X_best[4,:]

    KE_best = (1/2)*m_l*v_ls**2
    PE_best = (1/2)*k_test*(x_as*T1s - x_ls_adj/T2s)**2

    axs[0].plot(ts, KE_best,  linestyle='--', label=f"kinetic of best", alpha=0.8, linewidth=2)
    axs[0].plot(ts, PE_best,  linestyle='--', label=f"potential of best", alpha=0.8, linewidth=2)
    # ax0_twin.step(ts, X_best[4], 'b--', label=f"best", alpha=0.8, linewidth=2)
    # ax0_twin.set_ylim(0, 20)
    axs[0].set_xlabel('Time (seconds)')
    axs[0].set_ylabel('Energy (J)')
    # ax0_twin.set_ylabel('T2', color='blue')

    axs[0].legend(loc='upper left')
    # ax0_twin.legend(loc='lower right')

    plt.title(f"Kinetic and potential energy for top {topN} transmission profiles\nk={k_test}, m={m_l}, t={t_total}, L={leg_max}, v_best={round(v_max,3)}")

    os.makedirs(f"{exp_path}/figs_legvel", exist_ok=True)
    plt.savefig(f'{exp_path}/figs_legvel/{param_str}.png', dpi=300)
    plt.close()

    # break



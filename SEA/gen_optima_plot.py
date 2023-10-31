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
        
        if(v_file > v_max):
            v_max = v_file
            X_best = X_optload

        ts = np.linspace(0, t_total, nsteps+1)
        axs[0].step(ts, X_optload[2], 'r', alpha=0.2)
        # ax0_twin.step(ts, X_optload[4], 'b', alpha=0.2)

    axs[0].step(ts, X_best[2], 'r--', label=f"best", alpha=0.8, linewidth=2)
    # ax0_twin.step(ts, X_best[4], 'b--', label=f"best", alpha=0.8, linewidth=2)
    # ax0_twin.set_ylim(0, 20)
    axs[0].set_xlabel('Time (seconds)')
    axs[0].set_ylabel('Leg velociy (m/s)')
    # ax0_twin.set_ylabel('T2', color='blue')

    plt.title(f"Leg velocity for top {topN} transmission profiles\nk={k_test}, m={m_l}, t={t_total}, L={leg_max}, v_best={v_max}")

    os.makedirs(f"{exp_path}/figs_legvel", exist_ok=True)
    plt.savefig(f'{exp_path}/figs_legvel/{param_str}.png', dpi=300)
    plt.close()



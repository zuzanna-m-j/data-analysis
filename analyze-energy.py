#!/usr/bin/env python3

import glob
import os
import numpy as np
from celluloid import Camera
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.measure import label, regionprops, regionprops_table
from skimage.transform import rotate
import pickle
from matplotlib import cm
import matplotlib.colors as mcolors
import copy
from skimage.morphology import dilation,remove_small_holes, square, remove_small_objects, closing, disk
import matplotlib.colors
import sys
import matplotlib.ticker as mtick
import argparse
# plt.rcParams['figure.constrained_layout.use'] = True

parser = argparse.ArgumentParser()
parser.add_argument('--chi', default = "chips")
parser.add_argument('--dim', type=int, default=2)
parser.add_argument('--wlc', action='store_true')
parser.add_argument('--bmon', action='store_true')
args = parser.parse_args()


if args.chi == "chips":
    chi_name = r"$\chi_{PS}$"
elif args.chi == "chipi":
    chi_name = r"$\chi_{PI}$"
elif args.chi == "wlc":
    chi_name = r"$\wlc$"

colorsCS = list(mcolors.CSS4_COLORS)[10::5]
colors = ['orangered','tab:purple', 'dodgerblue', 'gold', 'forestgreen', 'darkred', 'darkorchid', 'darkorange', 'cornflowerblue', 'magenta', 'crimson', 'peru', 'olivedrab', 'gray']

rootdir = os.getcwd()
dirs = glob.glob(f'{rootdir}/[!_]*/')
dirs.sort()
with open(f"{rootdir}/data-summary.txt", "w") as f:
    f.writelines("chi_PS rho_p_rich rho_p_poor rho_ion_rich rho_ion_poor\n")

for dir in dirs:
    # file = glob.glob(f'{dir}data.dat')[0]
    file = f'{dir}data.dat'

    print(f"\n\nAnalyzing: {file}")
    if args.dim == 2:
        p_labels = ["xx", "yy"]
    else:
        p_labels = ["xx", "yy",'zz']

    file_name = "Blank" #file.split("/")[-3] + "-" + file.split("/")[-2]
    with open(f'{dir}input','r') as f_:
        lines = f_.readlines()
        for line in lines:
            items = line.strip().split()
            if len(items) > 0 and items[0] == 'n_gaussians':
                n_gauss = int(items[1])
                break

    print(n_gauss)
    data_file = np.loadtxt(file,skiprows = 3)

    col = 0

    time_step = data_file[:,col]
    col += 1

    # step Upe Ubond Angles Pxx Pyy Pxy
    UPe = data_file[:,col]
    col += 1

    UBond = data_file[:,col]
    col += 1

    if args.wlc:
        UAngle = data_file[:,col]
        col += 1

    Pressure = []
    if args.dim == 2:
        p_rng = 3
    else:
        p_rng = 6

    for i in range(col,col + p_rng):
        Pressure.append(data_file[:,i])

    col += p_rng

    UGauss = []
    for i in range(-n_gauss,0):
        UGauss.append(data_file[:,i])

    UGauss = np.array(UGauss)
    UBond = np.array(UBond)

    if args.wlc:
        UAngle = np.array(UAngle)

    UPe = np.array(UPe)
    UGaussSum = np.sum(UGauss, axis = 0)

    if args.wlc:
        UEl = UPe - UGaussSum -UBond - UAngle
    else:
        UEl = UPe - UGaussSum - UBond
        print("No wlc")
    if args.wlc:
        EN = [UPe, UEl, UAngle, UGaussSum, UBond]
        en_labels = ["UPe", "UEl", "UAngle", "UGaussSum", "UBon"]
    else:
        EN = [UPe, UEl, UGaussSum, UBond]
        en_labels = ["UPe", "UEl", "UGaussSum", "UBon"]

    fig, ax = plt.subplots(len(EN),1)
    for i, e in enumerate(EN):
        ax[i].plot(time_step, e, label = en_labels[i], color = colors[i])
        ax[i].legend()
        ax[i].set_yscale('log')
        if i != len(EN) - 1:
            ax[i].set_xticks([])
        box = ax[i].get_position()
        ax[i].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = 1, fancybox=True)
    ax[i].set_xlabel("Time Step")
    ax[i].set_ylabel("Energy")
    plt.tight_layout()
    


    # for i, p in enumerate(Pressure[:p_rng]):
    #     ax[1].plot(time_step,p, label = p_labels[i], color = colors[i])

    # ax[1].set_xlabel("Time Step")
    # ax[1].set_ylabel("Pressure")
    # ax[1].legend()

    plt.suptitle(f"{file_name}")
    plt.tight_layout()
    plt.show()
    #plt.savefig(f"histogram-clust.png", dpi = 300)

    fig, ax = plt.subplots(1,1)
    types = ["A","B","C","W","CAT","ANI","CI","D"]
    ax.plot(time_step,UGaussSum, label = 'UGaussSum', color = "dodgerblue")
    ax.set_xlabel("Time Step")
    ax.set_ylabel("Energy")
    ax.legend(loc = 'best')
    ax.set_yscale('log')
    plt.show()

    fig, ax = plt.subplots(1,1)

    g_labels = []
    for i in range(7):
        for j in range(i,7):
            g_labels.append(f"{types[i]} {types[j]}")
    c = 0
    for i,ugauss in enumerate(UGauss):
        if ugauss[0] != 0 and i not in (18,21):
            ax.plot(time_step,ugauss, label = g_labels[i], color = colors[c])
            c += 1
    ax.set_xlabel("Time Step")
    ax.set_ylabel("Energy")
    ax.legend(loc = 'best')
    ax.set_yscale('log')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol = 1, fancybox=True)
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
    #       fancybox=True, shadow=True, ncol=10)

    # plt.tight_layout()
    plt.show()

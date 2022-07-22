#!/usr/bin/env python3

from ast import arg
import glob
import os
from re import I
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

import numpy as np
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument('--chi', default = "chips")
parser.add_argument('--dim', default = 3, type=int)
parser.add_argument('-a', type = int, default=0)
parser.add_argument('-b', type = int, default=0)
parser.add_argument('-c', type = int, default=0)


args = parser.parse_args()

if args.chi == "chips":
    chi_name = r"$\chi_{PS}$"
elif args.chi == "chibs":
    chi_name = r"$\chi_{BS}$"
elif args.chi == "salt":
    chi_name = "salt ratio"
elif args.chi == "wlc":
    chi_name = "wlc"
elif args.chi == "bipolar":
    chi_name = r"$\alpha_{b}$"

dim = args.dim

# colors = list(mcolors.CSS4_COLORS)[10::3]
colors2 = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple','tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

#rootdir = '/home/jello/results/AB-block-21-Jun-2022/non-polar'

rootdir = os.getcwd()
polarity = os.path.basename(rootdir)
dirs = glob.glob(f'{rootdir}/[!_]*/')
dirs.sort()


CHI_PS = []
HIST_DATA = []
CLUST_SIZE = []
RHO_THR = []
DENS_PROF = [] 
ION_DENS_PROF = [] 
WATER_DENS_PROF = [] 

AC_CL_DENS = []
A_CL_DENS = []
B_CL_DENS = []

W_AC = []
W_A = []
W_B = []

I_AC = []
I_A = []
I_B = []


fig_anim, ax_anim = plt.subplots(1,1)
camera = Camera(fig_anim)


for dir in dirs[-1:]:

    REG_CL = []
    REG_CR = []

    hist_data = []

    chi_ps = float(dir[-5:-1])
    dir_name = f"{rootdir}"

    data = []

    files = glob.glob(f'{dir}*.tec')
    files.sort()

    files = files[-10:]
    file = files[0]
    print(f"\n\nAnalyzing: {file}")

    rho1_clusters = []
    rho1_clusters_ = []

    file_name = file.split("/")[-3] + "-" + file.split("/")[-2]

    f = np.loadtxt(file,skiprows = 3)
    xslices = list(set(f[:,0]))
    xslices.sort()

    yslices = list(set(f[:,1]))
    yslices.sort()
    yslices = np.array(yslices)

    if dim == 3:

        zslices = list(set(f[:,2]))
        zslices.sort()
        zslices = np.array(zslices)


######################################################################################################################

    if dim == 3:

        xdata = f[:,0]
        ydata = f[:,1]
        zdata = f[:,2]

        rhoA = f[:,3]/len(files)
        rhoB = f[:,4]/len(files)
        rhoC = f[:,5]/len(files)

        rhoW = f[:,6]/len(files)
        rhoCAT = f[:,7]/len(files)
        rhoANI = f[:,8]/len(files)
        rhoCI = f[:,9]/len(files)


        for data_file in files[1:]:

            f = np.loadtxt(data_file,skiprows = 3)

            rhoA += f[:,3]/len(files)
            rhoB += f[:,4]/len(files)
            rhoC += f[:,5]/len(files)

            rhoW += f[:,6]/len(files)
            rhoCAT += f[:,7]/len(files)
            rhoANI += f[:,8]/len(files)
            rhoCI += f[:,9]/len(files)

        data = [xdata,ydata,zdata,rhoA,rhoB,rhoC,rhoW,rhoCAT,rhoANI,rhoCI]


        b = 1.0
        N = args.a + args.b + args.c
        R_g = np.sqrt( (N-1) * (b**2)/6 )

        CHI_PS.append(chi_ps)

        rhoA_ = []
        rhoAC_ = []
        rhoB_ = []
        rhoIon_ = []
        rhoW_ = []

        x = xdata
        y = ydata
        z = zdata

        rhoAC = rhoA + rhoC
        rhoIon = rhoCAT + rhoANI + rhoCI

        for zslice in zslices:
            for i in range(len(z)):
                if z[i] == zslice:
                    rhoA_.append(rhoA[i])
                    rhoAC_.append(rhoAC[i])
                    rhoB_.append(rhoB[i])
                    rhoIon_.append(rhoIon[i])
                    rhoW_.append(rhoW[i])

        rhoAdata = np.array(rhoA_)
        rhoACdata = np.array(rhoAC_)
        rhoBdata = np.array(rhoB_)
        rhoIondata = np.array(rhoIon_)
        rhoWdata = np.array(rhoW_)

        rhoAdata *= (R_g**2/N)
        rhoACdata *= (R_g**2/N)
        rhoBdata *= (R_g**2/N)

        ion_norm = (args.a + args.c)/N * int(3 * 0.4)

        rhoIondata /= ion_norm

        ABC_dens = (rhoACdata + rhoBdata).reshape(len(zslices),-1)
        ABC_dens = np.mean(ABC_dens, axis = 1)
        DENS_PROF.append(ABC_dens)

        Ion_dens = rhoIondata.reshape(len(zslices),-1)
        Ion_dens = np.mean(Ion_dens, axis = 1)
        ION_DENS_PROF.append(Ion_dens)

        W_dens = rhoWdata.reshape(len(zslices),-1)
        W_dens = np.mean(W_dens, axis = 1)
        WATER_DENS_PROF.append(W_dens)



##############################################################################################

    elif dim == 2:

        xdata = f[:,0]
        ydata = f[:,1]

        rhoA = f[:,2]/len(files)
        rhoB = f[:,3]/len(files)
        rhoC = f[:,4]/len(files)

        rhoW = f[:,5]/len(files)
        rhoCAT = f[:,6]/len(files)
        rhoANI = f[:,7]/len(files)
        rhoCI = f[:,8]/len(files)


        for data_file in files[1:]:

            f = np.loadtxt(data_file,skiprows = 3)


            rhoA += f[:,2]/len(files)
            rhoB += f[:,3]/len(files)
            rhoC += f[:,4]/len(files)

            rhoW += f[:,5]/len(files)
            rhoCAT += f[:,6]/len(files)
            rhoANI += f[:,7]/len(files)
            rhoCI += f[:,8]/len(files)

        b = 1.0
        N = args.a + args.b + args.c
        R_g = np.sqrt( (N-1) * (b**2)/6 )

        CHI_PS.append(chi_ps)

        rhoA_ = []
        rhoAC_ = []
        rhoB_ = []
        rhoIon_ = []
        rhoW_ = []

        x = xdata
        y = ydata

        rhoAC = rhoA + rhoC
        rhoIon = rhoCAT + rhoANI + rhoCI

        for xslice in xslices:
            for i in range(len(x)):
                if x[i] == xslice:
                    rhoA_.append(rhoA[i])
                    rhoAC_.append(rhoAC[i])
                    rhoB_.append(rhoB[i])
                    rhoIon_.append(rhoIon[i])
                    rhoW_.append(rhoW[i])

        rhoAdata = np.array(rhoA_)
        rhoACdata = np.array(rhoAC_)
        rhoBdata = np.array(rhoB_)
        rhoIondata = np.array(rhoIon_)
        rhoWdata = np.array(rhoW_)


        rhoAdata *= (R_g**2/N)
        rhoACdata *= (R_g**2/N)
        rhoBdata *= (R_g**2/N)

        ion_norm = (args.a + args.c)/N * int(3 * 0.4)
        rhoIondata /= ion_norm

        #calculate density profile


        ABC_dens = (rhoACdata + rhoBdata).reshape(len(xslices),-1)
        ABC_dens = np.mean(ABC_dens, axis = 0)
        DENS_PROF.append(ABC_dens)

        Ion_dens = rhoIondata.reshape(len(xslices),-1)
        Ion_dens = np.mean(Ion_dens, axis = 0)
        ION_DENS_PROF.append(Ion_dens)

        W_dens = rhoWdata.reshape(len(xslices),-1)
        W_dens = np.mean(W_dens, axis = 0)
        WATER_DENS_PROF.append(W_dens)



#######################################################################################################################


    ac_dens = []
    a_dens = []
    b_dens = []

    w_ac = []
    w_a = []
    w_b = []

    i_ac = []
    i_a = []
    i_b = []


    A_thr = []
    for i in range(len(rhoACdata)):
        if rhoAdata[i] < 0.075:
            A_thr.append(0)
        else:
            A_thr.append(1)
            a_dens.append(rhoAdata[i])
            w_a.append(rhoWdata[i])
            i_a.append(rhoIondata[i])
    A_thr = np.array(A_thr)

    AC_thr = []
    for i in range(len(rhoACdata)):
        if rhoACdata[i] < 0.15:
            AC_thr.append(0)
        else:
            AC_thr.append(1)
            ac_dens.append(rhoACdata[i])
            w_ac.append(rhoWdata[i])
            i_ac.append(rhoIondata[i])
    AC_thr = np.array(AC_thr)

    B_thr = []
    for i in range(len(rhoBdata)):
        if rhoBdata[i] < 0.075:
            B_thr.append(0)
        else:
            B_thr.append(1)
            b_dens.append(rhoBdata[i])
            w_b.append(rhoWdata[i])
            i_b.append(rhoIondata[i])
    B_thr = np.array(B_thr)

    ABC_thr = np.zeros_like(B_thr)
    B_mask = np.argwhere(B_thr == 1)
    AC_mask = np.argwhere(AC_thr == 1)
    ABC_thr[B_mask] = 5
    ABC_thr[AC_mask] = 10


    AC_CL_DENS.append(np.average(ac_dens))
    A_CL_DENS.append(np.average(a_dens))
    B_CL_DENS.append(np.average(b_dens))

    W_AC.append(np.average(w_ac))
    W_A.append(np.average(w_a))
    W_B.append(np.average(w_b))

    I_AC.append(np.average(i_ac))
    I_A.append(np.average(i_a))
    I_B.append(np.average(i_b))

    if dim == 2:
        
        xdim, ydim = (AC_thr.reshape(len(xslices),-1)).shape #rows - x, cols - y

        footprint = square(12)
        # imgAC  = AC_thr.reshape(len(xslices),-1)
        # imgAC = closing(imgAC,footprint)
        # imgAC = remove_small_holes(imgAC, 12)
        # imgAC = remove_small_objects(imgAC, 5)
        # labels = label(imgAC)
        # regions = regionprops(labels)

        # imgB  = B_thr.reshape(len(xslices),-1)
        # imgB = closing(imgB,footprint)
        # imgB = remove_small_holes(imgB, 12)
        # imgB = remove_small_objects(imgB, 5)
        # labels = label(imgB)
        # regions = regionprops(labels)

        ax_anim.imshow(ABC_thr.reshape(len(xslices),-1).T, cmap = "viridis")
        ax_anim.set_xlim(0,xdim)
        ax_anim.set_ylim(0,ydim)
        ax_anim.text(200, 0.0, f"{chi_name}: {chi_ps}")
        ax_anim.set_xlabel('X')
        ax_anim.set_ylabel('Y')
        ax_anim.set_title("Polymer density")

        # ax_anim[0].imshow(ABC_thr.reshape(len(xslices),-1).T, cmap = "viridis")
        # ax_anim[0].set_xlim(0,xdim)
        # ax_anim[0].set_ylim(0,ydim)
        # ax_anim[0].text(200, 0.0, f"{chi_name}: {chi_ps}")
        # ax_anim[0].set_xlabel('X')
        # ax_anim[0].set_ylabel('Y')
        # ax_anim[0].set_title("Polymer density")

        # ax_anim[1].imshow(AC_thr.reshape(len(xslices),-1).T, cmap = "coolwarm")
        # ax_anim[1].set_xlim(0,xdim)
        # ax_anim[1].set_ylim(0,ydim)
        # ax_anim[1].set_xlabel('X')
        # ax_anim[1].set_ylabel('Y')
        # ax_anim[1].set_title("Threshold")

        fig_anim.suptitle(f"Polymer density ({polarity})")
        camera.snap()

        fig, ax = plt.subplots(1,3,constrained_layout = True)
        im0 = ax[0].imshow(rhoACdata.reshape(len(xslices),-1).T, cmap = "Reds")
        ax[0].set_xlim(0,xdim)
        ax[0].set_ylim(0,ydim)
        # ax[0].text(200, 0.0, f"{chi_name}: {chi_ps}")
        ax[0].set_xlabel('X')
        ax[0].set_ylabel('Y')
        ax[0].set_title("A & C")
        # fig.colorbar(im0, ax=ax[0])

        im1 = ax[1].imshow(ABC_thr.reshape(len(xslices),-1).T, cmap = "viridis")
        ax[1].set_xlim(0,xdim)
        ax[1].set_ylim(0,ydim)
        ax[1].set_xlabel('X')
        ax[1].set_ylabel('Y')
        # ax[1].set_title("Threshold")
        # fig.colorbar(im1, ax=ax[1])

        im1 = ax[2].imshow(B_thr.reshape(len(xslices),-1).T, cmap = "Greens")
        ax[2].set_xlim(0,xdim)
        ax[2].set_ylim(0,ydim)
        ax[2].set_xlabel('X')
        ax[2].set_ylabel('Y')
        ax[2].set_title("B")
        # fig.colorbar(im1, ax=ax[2])

        fig.suptitle(f"Polymer density - {polarity} - {chi_name}: {chi_ps}")
        fig.savefig(f"Polymer_density_{chi_ps}.png",dpi = 300)
        plt.close()

        # fig, ax = plt.subplots(1,3,constrained_layout = True)

        # im0 = ax[0].imshow(rhoIondata.reshape(len(xslices),-1).T, cmap = "Greens")
        # ax[0].set_xlim(0,xdim)
        # ax[0].set_ylim(0,ydim)
        # #ax[0].text(200, 0.0, f"{chi_name}: {chi_ps}")
        # ax[0].set_xlabel('X')
        # ax[0].set_ylabel('Y')
        # ax[0].set_title("Ions")
        # # fig.colorbar(im0, ax=ax[0])

        # im1 = ax[1].imshow(ABC_thr.reshape(len(xslices),-1).T, cmap = "viridis")
        # ax[1].set_xlim(0,xdim)
        # ax[1].set_ylim(0,ydim)
        # ax[1].set_xlabel('X')
        # ax[1].set_ylabel('Y')
        # ax[1].set_title("Polymer")
        # # fig.colorbar(im1, ax=ax[1])

        # im1 = ax[2].imshow(rhoWdata.reshape(len(xslices),-1).T, cmap = "Blues")
        # ax[2].set_xlim(0,xdim)
        # ax[2].set_ylim(0,ydim)
        # ax[2].set_xlabel('X')
        # ax[2].set_ylabel('Y')
        # ax[2].set_title("Water")
        # # fig.colorbar(im1, ax=ax[2])

        # fig.suptitle(f"Species density ({polarity})")
        # fig.savefig(f"Species_density_{chi_ps}.png",dpi = 500)
        # plt.close()


            # corona = ((dilation(img,square(7))).astype(int) - img.astype(int)).astype(bool)
            # corona_flat = ((dilation(img,square(7))).astype(int) - img.astype(int)).ravel()

        # Ion_thr = []

        # for i in range(len(rhoIondata)):
        #     if abs(rhoIondata[i] - meanIon) < stdIon:
        #         Ion_thr.append(0)
        #     else:
        #         Ion_thr.append(1)
        # Ion_thr = np.array(Ion_thr)

        # _p_AC_rich = []
        # _p_B_rich = []
        # _p_solv = []

        # _i_AC_rich = []
        # _i_B_rich = []
        # _i_solv = []



        # thr = np.argwhere(B_thr == 1)

        # for r1,r2 in zip(rhoACdata[thr],rhoIondata[thr]):
        #     _p_B_rich.append(r1)
        #     _i_B_rich.append(r2)


        # thr = np.argwhere((AC_thr + B_thr)== 0)
        # for r1,r2 in zip(rhoACdata[thr],rhoIondata[thr]):
        #     _p_solv.append(r1)
        #     _i_solv.append(r2)

        # i_AC_rich.append(np.average(_i_AC_rich))
        # i_B_rich.append(np.average(_i_B_rich))
        # i_solv.append(np.average(_i_solv))


if  dim == 3:
    yslices = zslices


max_rho = []
chi_rho = []

def fit1(x,cs,cp,squiggle,l0):
    return cs + (cp - cs) * 1/2 * (np.tanh((x-l0)/squiggle) + 1)

def fit2(x,cs,cp,squiggle,l0):
    return cs + (cp - cs) * 1/2 * (np.tanh((-x+l0)/squiggle) + 1)

fig, ax = plt.subplots()
plt.suptitle(f"Polymer Density Profile - Fitted - {polarity}")
ax.set_xlabel('y')
if dim == 3:
    ax.set_xlabel('z')
ax.set_ylabel(r"$<C^{*}>$")

for i in range(len(DENS_PROF)):
    try:
        y1 = yslices[:len(DENS_PROF[i])//2]
        vals1 = DENS_PROF[i][:len(DENS_PROF[i])//2]
        cs = vals1[20]
        cp = max(vals1)
        squiggle = 1.0
        l0 = 50
        popt1, pcov1 = curve_fit(fit1, y1,vals1, p0 = [cs,cp, squiggle, l0])

        y2 = yslices[len(DENS_PROF[i])//2:]
        vals2 = DENS_PROF[i][len(DENS_PROF[i])//2:]
        cs = vals2[-20]
        cp = max(vals2)
        squiggle = 1.0
        l0 = 80
        popt2, pcov2 = curve_fit(fit2, y2,vals2, p0 = [cs,cp, squiggle, l0])

        cp = (popt1[1] + popt2[1])/2

        # ax.plot(y1,vals1,marker = '.', linestyle = "None", alpha = 0.5, color = colors2[i])
        ax.plot(y1,fit1(y1,popt1[0], cp, *popt1[2:]),color = colors2[i], label = f'{chi_name}: {CHI_PS[i]}')

        # ax.plot(y2,vals2,marker = '.', linestyle = "None", alpha = 0.2, color = colors2[i])
        ax.plot(y2,fit2(y2,popt2[0], cp, *popt2[2:]),color = colors2[i])
        max_rho.append(max(fit1(y1,popt1[0], cp, *popt1[2:])))
        chi_rho.append(CHI_PS[i])

    except:
        print(f"Error: failed on {CHI_PS[i]}!\n")

with open('fit.data','w') as f:
    f.writelines(f"peak_rho = {str(max_rho)}\n")
    f.writelines(f"chi_vals = {str(chi_rho)}\n")


ax.legend(loc = 'best')
plt.tight_layout()
fig.savefig(f"FIT_POL_{polarity}.png", dpi = 300)
plt.close()


res = camera.animate(interval = 500)
res.save(f'density_profile_{polarity}.gif', dpi = 300)

CL_POL_DENS = []
CL_ION_DENS = []
CL_WATER_DENS = []

SOL_POL_DENS = []
SOL_ION_DENS = []



ordinate = yslices

fig, ax = plt.subplots()
plt.suptitle(f"Polymer Density Profile ({polarity})")
for i in range(len(DENS_PROF)):

    vals = DENS_PROF[i]
    ivals = ION_DENS_PROF[i]
    wvals = WATER_DENS_PROF[i]

    thr_high = np.argwhere(vals > 0.3)
    thr_low = np.argwhere(vals <= 0.3)

    CL_POL_DENS.append(np.average(vals[thr_high]))
    CL_ION_DENS.append(np.average(ivals[thr_high]))
    CL_WATER_DENS.append(np.average(wvals[thr_high]))

    SOL_POL_DENS.append(np.average(vals[thr_low]))
    SOL_ION_DENS.append(np.average(ivals[thr_low]))

    ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}')


ax.legend(loc = 'best')
ax.set_ylabel(r"$C^{*}$")
ax.set_xlabel("y")
plt.tight_layout()
fig.savefig(f"Polymer_Density_Profile_{polarity}.png", dpi = 300)


fig, ax = plt.subplots()
plt.suptitle(f"Polymer Density - {polarity}")
ax.set_ylabel(f"{chi_name}")
ax.set_xlabel(r"$<C^{*}>$")
ax.plot(CL_POL_DENS,CHI_PS, color = "tab:red", linestyle = 'None', marker = 'o', label = 'slab')
ax.plot(SOL_POL_DENS,CHI_PS, color = "tab:blue", linestyle = 'None', marker = 'd', label = 'solvent')
ax.legend(loc = 'best')
fig.savefig(f"Pol_Density_Clust_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
plt.suptitle(f"Ion Density Profile ({polarity})")
for i in range(len(ION_DENS_PROF)):
    vals = ION_DENS_PROF[i]
    ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
ax.legend(loc = 'best')
ax.set_ylabel(r"$\rho_{ion}$")
ax.set_xlabel("y")
plt.tight_layout()
fig.savefig(f"Ion_Density_Profile_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
plt.suptitle(f"Water Density Profile ({polarity})")
for i in range(len(WATER_DENS_PROF)):
    vals = WATER_DENS_PROF[i]
    ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
ax.legend(loc = 'best')
ax.set_ylabel(r"$\rho_{ion}$")
ax.set_xlabel("y")
plt.tight_layout()
fig.savefig(f"Water_Density_Profile_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
camera = Camera(fig)
plt.suptitle(f"Ion Density Profile ({polarity})")
for i in range(len(ION_DENS_PROF)):
    vals = ION_DENS_PROF[i]
    ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
    camera.snap()
ax.legend(loc = 'best')
ax.set_ylabel(r"$\rho_{ion}$")
ax.set_xlabel("y")
plt.tight_layout()
res = camera.animate(interval = 500)
res.save(f'Ion_Density_Profile_{polarity}.gif', dpi = 300)


fig, ax = plt.subplots()
plt.suptitle(f"Ion density - {polarity}")
ax.set_ylabel(f"{chi_name}")
ax.set_xlabel(r"$<\rho_{ion}>$")
ax.plot(CL_ION_DENS,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'o', label = 'slab')
ax.plot(SOL_ION_DENS,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'd', label = 'solvent')
ax.legend(loc = 'best')
fig.savefig(f"Ion_Density_Clust_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
plt.suptitle(f"Ion density - {polarity}")
ax.set_ylabel(f"{chi_name}")
ax.set_xlabel(r"$<\rho_{ion}>$")
ax.plot(CL_ION_DENS,CHI_PS, color = 'tab:green', linestyle = 'None', marker = 'o', label = 'slab')
ax.plot(SOL_ION_DENS,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'd', label = 'solvent')
ax.plot(I_AC,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'H', label = 'A&C')
# ax.plot(I_A,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'D', label = 'A')
ax.plot(I_B,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'X', label = 'B')
ax.legend(loc = 'best')
fig.savefig(f"Ion_Density_Clust_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
plt.suptitle(f"Monomer density - {polarity}")
ax.set_ylabel(f"{chi_name}")
ax.set_xlabel(r"$<C^{*}>$")
ax.plot(CL_POL_DENS,CHI_PS, color = "tab:green", linestyle = 'None', marker = 'o', label = 'slab')
ax.plot(A_CL_DENS,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'H', label = 'A')
ax.plot(AC_CL_DENS,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'D', label = 'A&C')
ax.plot(B_CL_DENS,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'X', label = 'B')
ax.legend(loc = 'best')
fig.savefig(f"Mon_Density_Clust_{polarity}.png", dpi = 300)
plt.close()

fig, ax = plt.subplots()
plt.suptitle(f"Water density - {polarity}")
ax.set_ylabel(f"{chi_name}")
ax.set_xlabel(r"$<C^{*}>$")
ax.plot(W_A,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'o', label = 'A')
ax.plot(W_AC,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'o', label = 'A&C')
ax.plot(W_B,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'd', label = 'B')
ax.legend(loc = 'best')
fig.savefig(f"Water_Density_Clust_{polarity}.png", dpi = 300)
plt.close()

with open('concentrations.data', 'w') as f:
    f.writelines(f'AC_conc = {str(AC_CL_DENS)}\n')
    f.writelines(f'A_conc = {str(A_CL_DENS)}\n')
    f.writelines(f'B_conc = {str(B_CL_DENS)}\n')
    f.writelines(f'IAC_conc = {str(I_AC)}\n')
    f.writelines(f'IA_conc = {str(I_A)}\n')
    f.writelines(f'IB_conc = {str(I_B)}\n')
    f.writelines(f'WAC_conc = {str(W_AC)}\n')
    f.writelines(f'WA_conc = {str(W_A)}\n')
    f.writelines(f'WB_conc = {str(W_B)}\n')




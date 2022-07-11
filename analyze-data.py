#!/usr/bin/env python3

from ast import arg
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

parser = argparse.ArgumentParser()
parser.add_argument('--chi', default = "chips")
parser.add_argument('--dim', default = 3, type=int)
parser.add_argument('-N', type = int)


args = parser.parse_args()

if args.chi == "chips":
    chi_name = r"$\chi_{PS}$"
elif args.chi == "chibs":
    chi_name = r"$\chi_{BS}$"
elif args.chi == "salt":
    chi_name = "salt ratio"

dim = args.dim

# colors = list(mcolors.CSS4_COLORS)[10::3]
colors2 = ['orangered','teal', 'dodgerblue', 'gold', 'forestgreen', 'darkred', 'darkorchid', 'darkorange', 'cornflowerblue', 'darkgreen', 'crimson', 'peru', 'olivedrab']

#rootdir = '/home/jello/results/AB-block-21-Jun-2022/non-polar'

rootdir = os.getcwd()
polarity = os.path.basename(rootdir)
dirs = glob.glob(f'{rootdir}/[!_]*/')
dirs.sort()

#
### 3D data arrays ###
#

polymer_rich = []
salt_rich = []
polymer_corona = []
salt_corona = []
polymer_poor = []
salt_poor = []

#
### 2D data arrays ###
#

CHI_PS = []
HIST_DATA = []
CLUST_SIZE = []
RHO_THR = []
DENS_PROF = [] 
ION_DENS_PROF = [] 
WATER_DENS_PROF = [] 

AC_CL_DENS = []
B_CL_DENS = []

W_AC = []
W_B = []

I_AC = []
I_B = []


fig_anim, ax_anim = plt.subplots(1,1)
camera = Camera(fig_anim)
       
with open(f"{rootdir}/summary.txt", "w") as f:
    f.writelines("chi_PS rho_p_rich rho_p_poor rho_ion_rich rho_ion_poor\n")

with open(f"{rootdir}/summary2.txt", "w") as f:
    f.writelines("chi_PS rho_pol_mean rho_salt_mean std_pol std_salt max_pol max_salt\n")

for dir in dirs:

    REG_CL = []
    REG_CR = []

    hist_data = []

    chi_ps = float(dir[-5:-1])
    dir_name = f"{rootdir}"

    #dir_name = f"{rootdir}/results"
    #os.system(f"mkdir {dir_name}")

    data = []

    files = glob.glob(f'{dir}*.tec')
    files.sort()

    file = files[-1]
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

    # load 2D data

    if dim == 2:

        xdata = f[:,0]
        ydata = f[:,1]

        rhoA = f[:,2]
        rhoB = f[:,3]
        rhoC = f[:,4]

        rhoW = f[:,5]
        rhoCAT = f[:,6]
        rhoANI = f[:,7]
        rhoCI = f[:,8]

        data = [xdata,ydata,rhoA,rhoB,rhoC,rhoW,rhoCAT,rhoANI,rhoCI]

    elif dim == 3:

        xdata = f[:,0]
        ydata = f[:,1]
        zdata = f[:,2]

        rhoA = f[:,3]
        rhoB = f[:,4]
        rhoC = f[:,5]

        rhoW = f[:,6]
        rhoCAT = f[:,7]
        rhoANI = f[:,8]
        rhoCI = f[:,9]

        data = [xdata,ydata,zdata,rhoA,rhoB,rhoC,rhoW,rhoCAT,rhoANI,rhoCI]

    # fig, ax = plt.subplots(1,4)
    # camera = Camera(fig)

    ydata = []
    zdata = []

    rho1data = []
    rho2data = []

    rhoACdata = []
    rhoBdata = []
    rhoIondata = []

    print("Collecting slice data...")

    # 2D data processing

    if dim == 2:

        N = args.N
        R_g = np.sqrt((N-1) * 1.0**2/6)

        CHI_PS.append(chi_ps)

        rhoAC_ = []
        rhoB_ = []
        rhoIon_ = []
        rhoW_ = []

        x = data[0]
        y = data[1]
        z = data[2]

        rhoAC = data[2] + data[4]
        rhoB = data[3]
        rhoW = data[5]
        rhoIon = data[6] + data[7] + data[8]

        for xslice in xslices:
            for i in range(len(x)):
                if x[i] == xslice:
                    rhoAC_.append(rhoAC[i])
                    rhoB_.append(rhoB[i])
                    rhoIon_.append(rhoIon[i])
                    rhoW_.append(rhoW[i])

        rhoACdata = np.array(rhoAC_)
        rhoBdata = np.array(rhoB_)
        rhoIondata = np.array(rhoIon_)
        rhoWdata = np.array(rhoW_)

        # meanAC = rhoACdata.mean()
        # meanB = rhoBdata.mean()
        # meanIon = rhoIondata.mean()

        rhoACdata *= (R_g**2/N)
        rhoBdata *= (R_g**2/N)
        # rhoIondata *= (R_g**2/N)

        meanAC = rhoACdata.mean()
        meanB = rhoBdata.mean()
        meanIon = rhoIondata.mean()

        stdAC = np.std(rhoACdata)
        stdB = np.std(rhoBdata)
        stdIon = np.std(rhoIondata)

        # print(f"Average polymer density: {mean1}, average CI density: {mean2}")
        # print(f"Max polymer density: {max1}, max CI density: {max2}")

        # with open(f"{rootdir}/summary2.txt", "a+") as f:
        #     f.writelines(f'{chi_ps} {mean1} {mean2} {std1} {std2} {max1} {max2}\n')  





        ### thresholding #######

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

        ac_dens = []
        b_dens = []

        w_ac = []
        w_b = []

        i_ac = []
        i_b = []

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
        B_CL_DENS.append(np.average(b_dens))

        W_AC.append(np.average(w_ac))
        W_B.append(np.average(w_b))

        I_AC.append(np.average(i_ac))
        I_B.append(np.average(i_b))

        
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

        fig, ax = plt.subplots(1,3,constrained_layout = True)

        im0 = ax[0].imshow(rhoIondata.reshape(len(xslices),-1).T, cmap = "Greens")
        ax[0].set_xlim(0,xdim)
        ax[0].set_ylim(0,ydim)
        #ax[0].text(200, 0.0, f"{chi_name}: {chi_ps}")
        ax[0].set_xlabel('X')
        ax[0].set_ylabel('Y')
        ax[0].set_title("Ions")
        # fig.colorbar(im0, ax=ax[0])

        im1 = ax[1].imshow(ABC_thr.reshape(len(xslices),-1).T, cmap = "viridis")
        ax[1].set_xlim(0,xdim)
        ax[1].set_ylim(0,ydim)
        ax[1].set_xlabel('X')
        ax[1].set_ylabel('Y')
        ax[1].set_title("Polymer")
        # fig.colorbar(im1, ax=ax[1])

        im1 = ax[2].imshow(rhoWdata.reshape(len(xslices),-1).T, cmap = "Blues")
        ax[2].set_xlim(0,xdim)
        ax[2].set_ylim(0,ydim)
        ax[2].set_xlabel('X')
        ax[2].set_ylabel('Y')
        ax[2].set_title("Water")
        # fig.colorbar(im1, ax=ax[2])

        fig.suptitle(f"Species density ({polarity})")
        fig.savefig(f"Species_density_{chi_ps}.png",dpi = 500)
        plt.close()


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


    elif dim == 3:
        print("You fucked up!")

        for xslice in xslices:

            rho1_ = []
            rho2_ = []

            x = data[0]
            y = data[1]
            z = data[2]

            rho1 = data[3] + data[4] + data[5] #polymer density
            rho2 = data[9] #counter-ion density

            for i in range(len(rho1)):
                if x[i] == xslice:
                    rho1_.append(rho1[i])
                    rho2_.append(rho2[i])
            rho1data.append(rho1_)
            rho2data.append(rho2_)


        rho1data = np.array(rho1data)
        rho2data = np.array(rho2data)

        # normalize the data

        mean1 = (rho1data.ravel()).mean()
        mean2 = (rho2data.ravel()).mean()

        rho1data /=mean1
        rho2data /=mean2

        m1 = (rho1data.ravel()).mean()
        m2 = (rho2data.ravel()).mean()

        max1 = (rho1data.ravel()).max()
        max2 = (rho2data.ravel()).max()

        std1 = np.std(rho1data.ravel())
        std2 = np.std(rho2data.ravel())

        print(f"Average polymer density: {mean1}, average CI density: {mean2}")
        print(f"Max polymer density: {max1}, max CI density: {max2}")

        with open(f"{rootdir}/summary2.txt", "a+") as f:
            f.writelines(f'{chi_ps} {mean1} {mean2} {std1} {std2} {max1} {max2}\n')


        rho1rich = []
        rho2rich = []

        rho2corona = []
        rho1corona = []

        rho1poor = []
        rho2poor = []

        #density correlation

        for x,xslice in enumerate(xslices):
            rho1_thr = []
            for i in range(len(rho1data[x])):
                if abs(rho1data[x][i] - m1) < std1:
                    rho1_thr.append(0)
                else:
                    rho1_thr.append(1)
            rho1_thr = np.array(rho1_thr)

            footprint = square(12)
            zdim, ydim = (rho1_thr.reshape(-1,len(xslices))).shape
            img  = rho1_thr.reshape(-1,len(xslices))
            img = closing(img,footprint)
            img = remove_small_holes(img, 12)
            img = remove_small_objects(img, 5)

            corona = ((dilation(img,square(7))).astype(int) - img.astype(int)).astype(bool)
            corona_flat = ((dilation(img,square(7))).astype(int) - img.astype(int)).ravel()
            labels = label(img)
            regions = regionprops(labels)

            rho2_thr = []

            for i in range(len(rho2data[x])):
                if abs(rho2data[x][i] - m2) < std2:
                    rho2_thr.append(0)
                else:
                    rho2_thr.append(1)
            rho2_thr = np.array(rho2_thr)

            img2 = rho2_thr.reshape(-1,len(xslices))

            img2_corona = img2.copy()
            img2_cluster = img2.copy()

            img2_corona[~corona] = 0
            img2_cluster[~img] = 0

            labels_cl = label(img2_cluster)
            regions_cl = regionprops(labels_cl)

            labels_cr = label(img2_corona)
            regions_cr = regionprops(labels_cr)

            REG_CL.append(regions_cl)
            REG_CR.append(regions_cr)

            dens_cl = np.zeros(len(rho1data[x]))
            for i in range(len(rho1data[x])):
                if rho1_thr[i] == 1:
                    if rho2_thr[i] == 1:
                        dens_cl[i] = 1
                    else:
                        dens_cl[i] = -1
                elif rho1_thr[i] == 0:
                    if rho2_thr[i] == 0:
                        dens_cl[i] = 1
                    else:
                        dens_cl[i] = -1


            dens_cr = np.zeros(len(rho2data[x]))
            for i in range(len(rho2data[x])):
                if corona_flat[i] == 1:
                    if rho2_thr[i] == 1:
                        dens_cr[i] = 1
                    else:
                        dens_cr[i] = -1
                elif corona_flat[i] == 0:
                    if rho2_thr[i] == 0:
                        dens_cr[i] = 1
                    else:
                        dens_cr[i] = -1


            for point in rho1data[x]:
                hist_data.append(point)
            
            # rho2corona = []
            # rho2cluster = []
            # rho1cluster = []

            # for t in range(len(rho1_thr)):
            #     if flat_corona[t] == 1:
            #         rho2corona.append([xslices[x], YDATA[x,t], ZDATA[x,t], rho1data[x,t]])
            #     else:
            #         rho2corona.append([xslices[x], YDATA[x,t], ZDATA[x,t], 0])
        
            #     if rho1_thr[t] == 1:
            #         rho1cluster.append([xslices[x], YDATA[x,t], ZDATA[x,t], rho1data[x,t]])
            #         rho2cluster.append([xslices[x], YDATA[x,t], ZDATA[x,t], rho1data[x,t]])
            #     else:
            #         rho1cluster.append([xslices[x], YDATA[x,t], ZDATA[x,t], 0])
            #         rho2cluster.append([xslices[x], YDATA[x,t], ZDATA[x,t], 0])

            # rho2corona = np.array(rho2corona)
            # rho2cluster = np.array(rho2cluster)
            # rho1cluster = np.array(rho1cluster)

            # Region preview - sanity check
            # fig, ax = plt.subplots()
            # ax.imshow(img, cmap=plt.cm.gray)

            # for j,props in enumerate(regions):

            #     minr, minc, maxr, maxc = props.bbox
            #     bx = (minc, maxc, maxc, minc, minc)
            #     by = (minr, minr, maxr, maxr, minr)
            #     ax.plot(bx, by, c = colors[j], linewidth=0.5)
            # ax.set_title("Selected clusters - 2D slice")
            # plt.show()

            # Diplay density data

            dens_cl = dens_cl.reshape(-1,len(xslices))
            dens_cr = dens_cr.reshape(-1,len(xslices))

            ax[0].imshow(img, cmap = "Greens")
            ax[0].set_xlim(0,ydim)
            ax[0].set_ylim(0,zdim)
            ax[0].set_title(r"Cluster")
            ax[0].text(ydim + 1, 0.0, f"X: {xslice}")
            ax[0].set_xlabel('Y')
            ax[0].set_ylabel('Z')

            ax[1].imshow(corona, cmap = "Blues")
            ax[1].set_xlim(0,ydim)
            ax[1].set_ylim(0,zdim)
            ax[1].set_title(r"Corona")
            #ax[1].text(ydim + 1, 0.0, f"X: {xslice}")
            # ax[1].set_xlabel('Y')
            # ax[1].set_ylabel('Z')

            ax[2].imshow(dens_cl, cmap = "coolwarm")
            ax[2].set_xlim(0,ydim)
            ax[2].set_ylim(0,zdim)
            ax[2].set_title(r"Corr. cluster")
            #ax[2].text(ydim + 1, 0.0, f"X: {xslice}")
            # ax[2].set_xlabel('Y')
            # ax[2].set_ylabel('Z')

            ax[3].imshow(dens_cr.reshape(-1,len(xslices)), cmap = "coolwarm")
            ax[3].set_xlim(0,ydim)
            ax[3].set_ylim(0,zdim)
            ax[3].set_title(r"Corr. corona")
            #ax[3].text(ydim + 1, 0.0, f"X: {xslice}")
            # ax[3].set_xlabel('Y')
            # ax[3].set_ylabel('Z')

            # ax[1].imshow(dens_cl.reshape(-1,len(xslices)), cmap = "Reds")
            # ax[1].set_xlim(0,ydim)
            # ax[1].set_ylim(0,zdim)
            # ax[1].set_title(r"Density correlation")
            # ax[1].text(ydim + 1, 0.0, f"X: {xslice}")
            # ax[1].set_xlabel('Y')
            # ax[1].set_ylabel('Z')

            # ax[3].imshow(dens_cl.reshape(-1,len(xslices)), cmap = "RdYlGn")
            # ax[3].set_xlim(0,ydim)
            # ax[3].set_ylim(0,zdim)
            # ax[3].set_title(r"Ion density")
            # ax[3].text(ydim + 1, 0.0, f"X: {xslice}")
            # ax[3].set_xlabel('Y')
            # ax[3].set_ylabel('Z')

            camera.snap()

            #ax.plot_trisurf(ydata,zdata,rho1data, cmap='viridis', edgecolor='none',alpha = 0.5)
            #ax.plot_surface(ydata[:-1].reshape(2,-1)[:-1].T, zdata[:-1].reshape(2,-1).T, rho1data[:-1].reshape(2,-1).T,cmap = 'viridis')
            #ax.scatter(ydata,zdata,rho1data, cmap='viridis', edgecolor='none',alpha = 0.5)


            # find region-specific concentration

            thr = np.argwhere(rho1_thr == 1)

            for r1,r2 in zip(rho1data[x][thr],rho2data[x][thr]):
                rho1rich.append(r1)
                rho2rich.append(r2)
            
            thr_c = np.argwhere(((corona.ravel()).astype(int)) == 1)
            for r1,r2 in zip(rho1data[x][thr_c], rho2data[x][thr_c]):
                rho1corona.append(r1)
                rho2corona.append(r2)

            thr2 = np.argwhere(rho1_thr == 0)
            for r1,r2 in zip(rho1data[x][thr2],rho2data[x][thr2]):
                rho1poor.append(r1)
                rho2poor.append(r2)

            if xslice == xslices[0]:
                print("Making clusters...")
                for region in regions:
                    rho1_clusters.append([[x,region.coords]])

            else:
                print(f"{x/len(xslices) * 100:3.1f}%", end = " ", flush=True)
                for regid,region in enumerate(regions):
                    overlap = False
                    neighs = []
                    coords = region.coords

                    for rid, p_region in enumerate(rho1_clusters):
                        if p_region != None:
                            for m in range(len(p_region)):
                                if overlap == True:
                                    break
                                dx = abs(x - p_region[m][0])
                                # if dx > len(xslices)/2:
                                #     dx = len(xslices) - dx
                                if dx == 1 or dx == 44:
                                    for coord in coords:
                                        if overlap == True:
                                            break
                                        dy  = abs(coord[1] - p_region[m][1][:,1])
                                        for i in range(len(dy)):
                                            if dy[i] >= int(ydim/2):
                                                dy[i] = ydim - dy[i]

                                        dz  = abs(coord[0] - p_region[m][1][:,0])
                                        for i in range(len(dz)):
                                            if dz[i] >= int(zdim/2):
                                                dz[i] = zdim - dz[i]
                                        dr = dy + dz
                                        if min(dr) <= np.sqrt(2):
                                            overlap = True
                                if overlap == True:
                                    neighs.append(rid)
                                    overlap = False
                                else:
                                    pass
                        
                    if len(neighs) != 0:
                        i = neighs[0]
                        for j in neighs:
                            if i != j:
                                cluster = rho1_clusters[j]
                                for slice in cluster:
                                    rho1_clusters[i].append(slice)

                        rho1_clusters[i].append([x, region.coords])
                        for j in neighs:
                            if i != j:
                                rho1_clusters[j] = None
                    else:
                            rho1_clusters.append([[x, region.coords]])

        #     RHO2CR.append(rho2corona)
        #     RHO1CL.append(rho1cluster)
        #     RHO2CL.append(rho1cluster)

        # RHO2CR = np.array(RHO2CR)
        # RHO1CL = np.array(RHO1CL)
        # RHO2CL = np.array(RHO2CL)

        HIST_DATA.append(hist_data)
        CHI_PS.append(chi_ps)

        print("\n")

        res = camera.animate(interval = 500)
        res.save(f'{dir_name}/rho{str(chi_ps)}.gif', dpi = 300)

        # plt.show()

        rho1rich = np.array(rho1rich)
        rho1rich_avg  = np.average(rho1rich)

        rho2rich = np.array(rho2rich)
        rho2rich_avg  = np.average(rho2rich)

        rho2corona = np.array(rho2corona)
        rho2corona_avg = np.average(rho2corona)

        rho1corona = np.array(rho1corona)
        rho1corona_avg = np.average(rho1corona)

        rho1poor = np.array(rho1poor)
        rho1poor_avg  = np.average(rho1poor)

        rho2poor = np.array(rho2poor)
        rho2poor_avg  = np.average(rho2poor)

        with open(f"{rootdir}/summary.txt", "a+") as f:
            f.writelines(f'{chi_ps} {(rho1rich_avg-m1)/std1} {(rho1poor_avg-m1)/std1} {(rho2rich_avg - m2)/std2} {(rho2poor_avg - m2)/std2}\n')

        polymer_rich.append((rho1rich_avg-m1)/std1)

        polymer_corona.append((rho1corona_avg - m1)/std1)

        salt_rich.append((rho2rich_avg - m2)/std2)

        salt_corona.append((rho2corona_avg - m2)/std2)

        polymer_poor.append((rho1poor_avg-m1)/std1)
        salt_poor.append((rho2poor_avg - m2)/std2)
        
        rho_file = open(f"{dir_name}/rho1_clusters2.pkl", "wb")
        pickle.dump(rho1_clusters,rho_file)
        rho_file.close()

        rho_file = open(f"{dir_name}/rho1_clusters2.pkl", "rb")
        clusters = pickle.load(rho_file)
        rho_file.close()

        
        xcl = []
        ycl = []
        zcl = []
        for m,slice in enumerate(REG_CL):
            for region in slice:
                for coord in region.coords:
                    xcl.append(m)
                    ycl.append(coord[1])
                    zcl.append(coord[0])
        voxel_cl = np.zeros((len(xslices),ydim,zdim))
        for i in range(len(xcl)):
            voxel_cl[xcl[i],ycl[i],zcl[i]] = 1

        xcr = []
        ycr = []
        zcr = []
        for m,slice in enumerate(REG_CR):
            for region in slice:
                for coord in region.coords:
                    xcr.append(m)
                    ycr.append(coord[1])
                    zcr.append(coord[0])
        voxel_cr = np.zeros((len(xslices),ydim,zdim))
        for i in range(len(xcr)):
            voxel_cr[xcr[i],ycr[i],zcr[i]] = 1


        colors  = ['orangered','teal', 'dodgerblue', 'gold', 'forestgreen', 'darkred', 'darkorchid', 'darkorange', 'cornflowerblue', 'darkgreen', 'crimson', 'peru', 'olivedrab']
        alpha = 1.0
        colors_rgb = []
        for color in colors:
            c = matplotlib.colors.to_rgba(color)
            colors_rgb.append((c[0],c[1],c[2],alpha))
        colors = colors_rgb

        cl_list = np.zeros(len(clusters))
        VOX = []

        print("Rendering clusters...")
        j = 0
        cluster_sizes = []
        fig = plt.figure(dpi = 1000)
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        for c,cluster in enumerate(clusters):
            if cluster != None:
                x = []
                y = []
                z = []

                for slice in cluster:

                    x_ = slice[0]
                    y_ = slice[1][:,1]
                    z_ = slice[1][:,0]
                    for i in range(len(y_)):
                        x.append(x_)
                        y.append(y_[i])
                        z.append(z_[i])

                if len(z)>500:
                    
                    cluster_sizes.append(len(z))

                    voxel = np.zeros((len(xslices),ydim,zdim))
                    for i in range(len(x)):
                        voxel[x[i],y[i],z[i]] = 1
                    VOX.append(copy.deepcopy(voxel))
                    color = colors[j]
                    ax.voxels(voxel, facecolors = color)
                    # ax.scatter(x,y,z,s = 0.2, c = colors[j], alpha = 1.0, label = f"{j}")
                    #ax.scatter(xi,zi,yi,s = 10.0, c = 'forestgreen', alpha = 0.25, label = f"{j}", norm = plt.Normalize(0, 10))

                    ax.set_xlabel("X")
                    ax.set_ylabel("Y")
                    ax.set_zlabel("Z")
                    ax.set_box_aspect((len(xslices),ydim,zdim))  # aspect ratio is 1:1:1 in data space
                    j += 1

        plt.title("Polymer density voxels")
        plt.savefig(f"{dir_name}/voxels-{str(chi_ps)}.png", dpi = 500)
        CLUST_SIZE.append(cluster_sizes)


        colors  = ['orangered','teal', 'dodgerblue', 'gold', 'forestgreen', 'darkred', 'darkorchid', 'darkorange', 'cornflowerblue', 'darkgreen', 'crimson', 'peru', 'olivedrab']
        alpha = 0.1
        colors_rgb = []
        for color in colors:
            c = matplotlib.colors.to_rgba(color)
            colors_rgb.append((c[0],c[1],c[2],alpha))
        colors = colors_rgb

        print("Rendering clusters with ions...")
        fig = plt.figure(dpi = 1000)
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        for vox in VOX:
            ax.voxels(vox, facecolors = colors[3])
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_box_aspect((len(xslices),ydim,zdim))  # aspect ratio is 1:1:1 in data space

        ax.voxels(voxel_cl, facecolors =  'cyan')
        ax.voxels(voxel_cr, facecolors =  'magenta')

        plt.title("Ion density voxels")
        plt.savefig(f"{dir_name}/voxels-ions-{str(chi_ps)}.png", dpi = 1000)

        print("Concentration plots...")
        fig, ax = plt.subplots()
        plt.suptitle("Density distribution")

        ax.set_title("Polymer concentration")
        ax.plot(CHI_PS, polymer_rich, 'o-', label = "cluster", color = 'tab:green')
        ax.plot(CHI_PS, polymer_poor, 'o-', label = "solvent", color = 'tab:blue')
        ax.plot(CHI_PS, polymer_corona, 'o-', label = "corona", color = 'tab:red')
        ax.set_xlabel(chi_name)
        ax.set_ylabel(r"$\frac{\rho_{pol} - <\rho_{pol}>}{\sigma_{pol}}$")
        ax.legend(loc = 'best')

        plt.tight_layout()
        plt.savefig(fname = f"{dir_name}/density-pol.png",dpi = 300)

        ## end 3d analysis


if dim == 3:

    fig, ax = plt.subplots()
    plt.suptitle("Density distribution")

    ax.set_title("Salt concentration")
    ax.plot(CHI_PS, salt_rich, 'o-', label = "cluster", color = 'tab:green')
    ax.plot(CHI_PS, salt_poor, 'o-', label = "solvent", color = 'tab:blue')
    ax.plot(CHI_PS, salt_corona, 'o-', label = "corona", color = 'tab:red')
    ax.plot(CHI_PS, np.array(salt_corona) + np.array(salt_rich), 'o-', label = "cl. + cr.", color = 'tab:purple')
    ax.set_xlabel(chi_name)
    ax.set_ylabel(r"$\frac{\rho_{salt} - <\rho_{salt}>}{\sigma_{salt}}$")
    ax.legend(loc = 'best')

    plt.tight_layout()
    plt.savefig(fname = f"{dir_name}/density-ion.png",dpi = 300)

    weights = []
    fig, ax = plt.subplots()
    for i in range(len(HIST_DATA)):
        weights.append(np.ones(len(HIST_DATA[i]))/len(HIST_DATA[i]))
    ax.hist(HIST_DATA, bins = 30, color=colors2[:len(HIST_DATA)], label=CHI_PS, weights=weights)
    ax.set_xlabel(r"$\frac{\rho}{\rho_{mean}}$")
    ax.set_ylabel("Volume fraction")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
    ax.set_yscale('log')
    plt.suptitle("Polymer clusters")
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"histogram.png", dpi = 300)

    h_file = open(f"{dir_name}/hist.pkl", "wb")
    pickle.dump(HIST_DATA,h_file)
    h_file.close()

    fig, ax = plt.subplots()
    ax.hist(CLUST_SIZE, bins = 20, color=colors2[:len(CLUST_SIZE)], label=CHI_PS)
    ax.set_xlabel(r"$\frac{\rho}{\rho_{mean}}$")
    ax.set_ylabel("Cluster size")
    ax.set_yscale('log')
    plt.suptitle("Polymer clusters")
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"histogram-clust.png", dpi = 300)

    h_file = open(f"{dir_name}/hist-clust.pkl", "wb")
    pickle.dump(HIST_DATA,h_file)
    h_file.close()

# fig, ax = plt.subplots()

# plt.suptitle("Cluster size distribution")

# ax[0].set_title("Cluster sizes")
# ax[0].plot(CHI_PS, polymer_rich, 'o-', label = "polymer-rich phase")
# ax[0].plot(CHI_PS, polymer_poor, 'o-', label = "polymer-poor phase")
# ax[0].xlabel(r"$\chi_{PS}$")
# ax[0].ylabel(r"$\frac{\rho_{pol} - <\rho_{pol}>}{\sigma_{pol}}$")

# ax[1].set_title("Average cluster size")
# ax[1].plot(CHI_PS, salt_rich, 'o-', label = "Salt in polymer-rich phase")
# ax[1].plot(CHI_PS, salt_rich, 'o-', label = "Salt in polymer-rich phase")
# ax[1].xlabel(r"$\chi_{PS}$")
# ax[1].ylabel(r"$\frac{\rho_{salt} - <\rho_{salt}>}{\sigma_{salt}}$")

# plt.tight_layout()
# plt.legend()
# plt.savefig(fname = f"{dir_name}/density.png",dpi = 1000)


# data analysis - see if the ions partition

# 1. How to detect polymer-rich phase
# 2. Characterize the density
# 3. Data correlation 0.05 0.05 0.20

# Plot points in 3D C = rho


# polymer rich phase - x - slice
#
# For each slice get average polymer density


    # A = [1,+1.0, 1] #charged +
    # B = [2, 0, 0] #neutral
    # C = [3,-1.0, 1] #charged -
    # W = [4, 0, 1] #non-polar
    # CAT = [5,+1,0]
    # ANI = [6,-1,0]
    # CI = [7,1,0]
    # D = [8, 1, 0] #drude oscillator


# analyze energies, etc - pull the data file
# plot
# do some data analysis on that --

# spatial correlations, deviation from the mean - pointwise - how - 


# data_file = np.loadtxt('/home/jello/results/AB-block-21-Jun-2022-old/non-polar/np1.00/data.dat',skiprows = 2)
# time_step = data_file[:,0]
# UPe = data_file[:,1]
# UBond = data_file[:,2]

# Pressure = []
# for i in range(3,9):
#     Pressure.append(data_file[:,i])

# UGauss = []
# for i in range(9, 9+28):
#     UGauss.append(data_file[:,i])

# # plt.plot(time_step, UPe, label = 'UPe')
# # plt.legend()
# # #plt.show()
# # plt.plot(time_step, UBond, label = 'Ubond')
# # plt.legend()
# # #plt.show()

# types = ["A","B","C","W","CAT","ANI","CI","D"]

# count = 0
# for i in range(7):
#     for j in range(i,7):
#         print(f"{count}: {types[i]} {types[j]}")
#         count += 1

# for i,ugass in enumerate(UGauss):
#     if ugass[0] != 0 and i not in (18,21):
#         plt.plot(time_step,ugass, label = f"{i}")
# plt.legend()
# #plt.show()

# for i,p in enumerate(Pressure[:3]):
#     plt.plot(time_step,p, label = f"{i}")
# plt.legend()
# #plt.show()

# 18 is water-water -- solvent!

# they seem to condense because it minimizes the energy of 

# detect clusters
# measure cluser size



# to do -- cluster size detection


# cluster size vs cluser shape
# electrostatic energy ---





# Analysis to do:

# make a dedicated folder
# make an array of data
# chi p_rich p_poor s_rich s_poor
# plot density fileds per slice and the - go forward and backwards
# plot polymer concentrations and salt concentration as we vary chi - save the plots

# analyze pressure and energy over the course of simulation - extract the electrostatic energy
# cluster size distribution histogram
# density distribution histogram





if dim == 2:

    # i_AC_rich = np.array(i_AC_rich)
    # i_B_rich = np.array(i_B_rich)
    # i_solv = np.array(i_solv)
    # i_POL = i_AC_rich + i_B_rich
    
    # fig, ax = plt.subplots()
    # plt.suptitle("Density distribution")

    # ax.set_title("Salt concentration")
    # ax.plot(CHI_PS, i_AC_rich, 'o-', label = "AC", color = 'tab:green')
    # ax.plot(CHI_PS, i_B_rich, 'o-', label = "B", color = 'tab:blue')
    # ax.plot(CHI_PS,i_POL, 'o-', label = "ABC", color = 'tab:red')
    # ax.plot(CHI_PS,i_solv, 'o-', label = "solvent", color = 'tab:purple')
    # ax.set_xlabel(chi_name)
    # ax.set_ylabel(r"$\frac{\rho_{salt} - <\rho_{salt}>}{\sigma_{salt}}$")
    # ax.legend(loc = 'best')

    # plt.tight_layout()
    # plt.show()


    # process data

    res = camera.animate(interval = 500)
    res.save(f'density_profile_{polarity}.gif', dpi = 300)

    CL_POL_DENS = []
    CL_ION_DENS = []
    CL_WATER_DENS = []

    SOL_POL_DENS = []
    SOL_ION_DENS = []

    ordinate = np.mean(yslices.reshape(-1, 11), axis=1)
    fig, ax = plt.subplots()
    plt.suptitle(f"Polymer Density Profile ({polarity})")
    for i in range(len(DENS_PROF)):

        vals = np.mean(DENS_PROF[i].reshape(-1, 11), axis=1)
        ivals = np.mean(ION_DENS_PROF[i].reshape(-1, 11), axis=1)
        wvals = np.mean(WATER_DENS_PROF[i].reshape(-1, 11), axis=1)
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

    ordinate = np.mean(yslices.reshape(-1, 11), axis=1)
    fig, ax = plt.subplots()
    plt.suptitle(f"Ion Density Profile ({polarity})")
    for i in range(len(ION_DENS_PROF)):
        vals = np.mean(ION_DENS_PROF[i].reshape(-1, 11), axis=1)
        ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
    ax.legend(loc = 'best')
    ax.set_ylabel(r"$\rho_{ion}$")
    ax.set_xlabel("y")
    plt.tight_layout()
    fig.savefig(f"Ion_Density_Profile_{polarity}.png", dpi = 300)
    plt.close()

    ordinate = np.mean(yslices.reshape(-1, 11), axis=1)
    fig, ax = plt.subplots()
    plt.suptitle(f"Ion Density Profile ({polarity})")
    for i in range(len(WATER_DENS_PROF)):
        vals = np.mean(WATER_DENS_PROF[i].reshape(-1, 11), axis=1)
        ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
    ax.legend(loc = 'best')
    ax.set_ylabel(r"$\rho_{ion}$")
    ax.set_xlabel("y")
    plt.tight_layout()
    fig.savefig(f"Water_Density_Profile_{polarity}.png", dpi = 300)
    plt.close()

    ordinate = np.mean(yslices.reshape(-1, 11), axis=1)
    fig, ax = plt.subplots()
    camera = Camera(fig)
    plt.suptitle(f"Ion Density Profile ({polarity})")
    for i in range(len(ION_DENS_PROF)):
        vals = np.mean(ION_DENS_PROF[i].reshape(-1, 11), axis=1)
        ax.plot(ordinate,vals,label = f'{chi_name}: {CHI_PS[i]}', marker = '.', linewidth = 1.0)
        camera.snap()
    ax.legend(loc = 'best')
    ax.set_ylabel(r"$\rho_{ion}$")
    ax.set_xlabel("y")
    plt.tight_layout()
    res = camera.animate(interval = 500)
    res.save(f'Ion_Density_Profile_{polarity}.gif', dpi = 300)


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
    ax.plot(CL_ION_DENS,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'o', label = 'slab')
    ax.plot(SOL_ION_DENS,CHI_PS, color = 'tab:blue', linestyle = 'None', marker = 'd', label = 'solvent')
    ax.plot(I_AC,CHI_PS, color = 'tab:green', linestyle = 'None', marker = 'H', label = 'A&C')
    ax.plot(I_B,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'X', label = 'B')
    ax.legend(loc = 'best')
    fig.savefig(f"Ion_Density_Clust_{polarity}.png", dpi = 300)
    plt.close()

    fig, ax = plt.subplots()
    plt.suptitle(f"Monomer density - {polarity}")
    ax.set_ylabel(f"{chi_name}")
    ax.set_xlabel(r"$<C^{*}>$")
    ax.plot(AC_CL_DENS,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'o', label = 'A&C')
    ax.plot(B_CL_DENS,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'd', label = 'B')
    ax.legend(loc = 'best')
    fig.savefig(f"Mon_Density_Clust_{polarity}.png", dpi = 300)
    plt.close()

    fig, ax = plt.subplots()
    plt.suptitle(f"Water density - {polarity}")
    ax.set_ylabel(f"{chi_name}")
    ax.set_xlabel(r"$<C^{*}>$")
    ax.plot(W_AC,CHI_PS, color = 'tab:red', linestyle = 'None', marker = 'o', label = 'A&C')
    ax.plot(W_B,CHI_PS, color = 'tab:purple', linestyle = 'None', marker = 'd', label = 'B')
    ax.legend(loc = 'best')
    fig.savefig(f"Water_Density_Clust_{polarity}.png", dpi = 300)
    plt.close()
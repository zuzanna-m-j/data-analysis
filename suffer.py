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








###############################!!!!!!!!!!!!!!!!!!!!!!!!#########################
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

fig, ax = plt.subplots()

plt.suptitle("Cluster size distribution")

ax[0].set_title("Cluster sizes")
ax[0].plot(CHI_PS, polymer_rich, 'o-', label = "polymer-rich phase")
ax[0].plot(CHI_PS, polymer_poor, 'o-', label = "polymer-poor phase")
ax[0].xlabel(r"$\chi_{PS}$")
ax[0].ylabel(r"$\frac{\rho_{pol} - <\rho_{pol}>}{\sigma_{pol}}$")

ax[1].set_title("Average cluster size")
ax[1].plot(CHI_PS, salt_rich, 'o-', label = "Salt in polymer-rich phase")
ax[1].plot(CHI_PS, salt_rich, 'o-', label = "Salt in polymer-rich phase")
ax[1].xlabel(r"$\chi_{PS}$")
ax[1].ylabel(r"$\frac{\rho_{salt} - <\rho_{salt}>}{\sigma_{salt}}$")

plt.tight_layout()
plt.legend()
plt.savefig(fname = f"{dir_name}/density.png",dpi = 1000)
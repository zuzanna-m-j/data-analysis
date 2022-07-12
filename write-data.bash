#!/bin/bash

#for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
#for salt in 0.005 0.01 0.05 0.10 0.15 0.2 0.25 0.30; do

clear
echo "Generating input files..."

cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/spice-up.py spice.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/write-data.bash write.bash

mkdir polar
cd polar

  python3 ../input.py --chips 1.0 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --polar --scramble --max_step 7500001 --traj_freq 500000 --midpush 0.05
  python3 ../stitch.py

  ~/Downloads/cuda-tild/gpu-tild

  for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
  mkdir "$chi"
  cd "$chi"
  python3 ../../input.py --chips "$chi" --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --polar
  python3 ../../stitch.py --atoms ../traj.lammpstrj
    cd ../
    done
cd ../


# cd polar
#   for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
#   cd "$chi"
#    ~/Downloads/cuda-tild/gpu-tild
#     cd ../
#     done
# cd ../



# clear
# echo "Generating input files..."


# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/spice-up.py spice.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/write-data.bash write.bash

# mkdir polar
# cd polar
#   for salt in 0.005 0.01 0.05 0.10 0.15 0.2 0.25 0.30; do
#   cd "$salt"
#   python3 ../../input.py --chips 0.75 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --salt "$salt" --polar
#   python3 ../../stitch.py --atoms traj.lammpstrj
#     cd ../
#     done
# cd ../


# mkdir non-polar
# cd non-polar

#   python3 ../input.py --chips 1.50 --chibs 2.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.05 --max_step 750001 --log_freq 5000 --traj_freq 50000 --phi 0.4
#   python3 ../stitch.py
#   ~/Downloads/cuda-tild/gpu-tild

#   for salt in 0.00 0.005 0.01 0.05 0.1 0.15 0.2 0.25 0.30 0.40; do
#   mkdir "$salt"
#   cd "$salt"
#   python3 ../../input.py --chips 1.00 --chibs 2.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --salt "$salt"
#   python3 ../../spice.py
#       cd ../
#     done
# cd ../





# rm -r AB-3d-slab-no-salt-Jul-10
# mkdir AB-3d-slab-no-salt-Jul-10
# cd AB-3d-slab-no-salt-Jul-10

# touch info.txt
# echo $'AC 3d system, this one is for equilibriation\n' >> info.txt
# echo $'No salt added to the system\n' >> info.txt

# #echo $'Folder polar: vary the chi_ps parameter\n' >> info.txt
# #echo $'Folder non-polar: vary the chi_ps parameter\n' >> info.txt
# #echo $'Folder counter-ions: chi_ps = 1.0, chi_pi varied\n' >> info.txt

# echo $'\n\n' >> info.txt

# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/write-data.bash .

# mkdir polar
# cd polar

#   ../input.py --chips 2.00 --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --midpush 0.05 --polar --max_step 70001 --log_freq 2000 --traj_freq 5000 --phi 0.3
#   ../stitch.py

  #~/Downloads/cuda-tild/gpu-tild

  # for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
  # mkdir "$chi"
  # cd "$chi"
  # ../../input.py --chips "$chi" --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --polar --phi 0.3
  # ../../stitch.py --atoms ../traj.lammpstrj
  #~/Downloads/cuda-tild/gpu-tild
#   cd ../
#   done
# cd ../

# mkdir non-polar
# cd non-polar

#   ../input.py --chips 2.00 --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --midpush 0.05 --max_step 70001 --log_freq 2000 --traj_freq 5000 --phi 0.3
#   ../stitch.py

  #~/Downloads/cuda-tild/gpu-tild

  # for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
  # mkdir "$chi"
  # cd "$chi"
  # ../../input.py --chips "$chi" --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --phi 0.3
  # ../../stitch.py --atoms ../traj.lammpstrj
  #~/Downloads/cuda-tild/gpu-tild
  # cd ../
  # done
cd ../


#
#------------------------------------- 2D - AB -------------------------------------------------
#


# rm -r AB-2d-salt-Jul-10-eq
# mkdir AB-2d-no-salt-Jul-9-eq
# cd AB-2d-no-salt-Jul-9-wq

# touch info.txt
# echo $'AB 2d system, test the role of the parameter chi\n' >> info.txt
# echo $'with salt\n' >> info.txt

#echo $'Folder polar: vary the chi_ps parameter\n' >> info.txt
#echo $'Folder non-polar: vary the chi_ps parameter\n' >> info.txt
#echo $'Folder counter-ions: chi_ps = 1.0, chi_pi varied\n' >> info.txt

# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py

# echo $'\n\n' >> info.txt

# mkdir polar
# cd polar

#   ../input.py --chips 0.75 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.05 --polar --max_step 10001 --log_freq 1000 --traj_freq 5000 --phi 0.4 --skip_solv
#   mkdir info
#   mv info.data info/
#   ../stitch.py
#   ~/Downloads/cuda-tild/gpu-tild
#   ../input.py --chips 0.75 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.03 --polar --max_step 15001 --log_freq 1000 --traj_freq 5000 --phi 0.4
#   ../stitch-eq.py --atoms traj.lammpstrj --info info/info.data
#   ~/Downloads/cuda-tild/gpu-tild


#   for salt in 0.00 0.005 0.01 0.05 0.1 0.15 0.2 0.25 0.30 0.40; do
#   mkdir "$salt"
#   cd "$salt"
#   ../../input.py --chips 1.5 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --polar --salt "$salt"
#   ../../stitch-eq.py --atoms ../traj.lammpstrj --info ../info.data
#   cd ../
#   done
# cd ../

# mkdir non-polar
# cd non-polar

#   ../input.py --chips 0.75 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.05 --max_step 10001 --log_freq 1000 --traj_freq 5000 --phi 0.4 --skip_solv
#   mkdir info
#   mv info.data info/
#   ../stitch.py
#   ~/Downloads/cuda-tild/gpu-tild
#   ../input.py --chips 0.75 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.03 --max_step 15001 --log_freq 1000 --traj_freq 5000 --phi 0.4
#   ../stitch-eq.py --atoms traj.lammpstrj --info info/info.data
#   ~/Downloads/cuda-tild/gpu-tild


#   for salt in 0.00 0.005 0.01 0.05 0.1 0.15 0.2 0.25 0.30 0.40; do
#   mkdir "$salt"
#   cd "$salt"
#   ../../input.py --chips 1.5 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --salt "$salt"
#   ../../stitch-eq.py --atoms ../traj.lammpstrj --info ../info.data
#   cd ../
#   done
# cd ../



#
#------------------------------------- 2D - ABC -------------------------------------------------
#


# rm -r ABC-2d-no-salt-Jul-10
# mkdir ABC-2d-no-salt-Jul-10
# cd ABC-2d-no-salt-Jul-10

# touch info.txt
# echo $'ABC 2d system, test the role of the parameter chi\n' >> info.txt
# echo $'No salt added to the system\n' >> info.txt

# # #echo $'Folder polar: vary the chi_ps parameter\n' >> info.txt
# # #echo $'Folder non-polar: vary the chi_ps parameter\n' >> info.txt
# # #echo $'Folder counter-ions: chi_ps = 1.0, chi_pi varied\n' >> info.txt

# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
# cp /home/jello/results/submits/write-data.bash .

# echo $'\n\n' >> info.txt

# mkdir polar
# cd polar

#   ../input.py --chips 1.00 --chibs 2.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.05 --polar --max_step 450001 --log_freq 1000 --traj_freq 10000 --phi 0.4 --bmon
#   ../stitch.py
#   # ../stitch.py --atoms traj.lammpstrj
#   ~/Downloads/cuda-tild/gpu-tild

#   for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
#   mkdir "$chi"
#   cd "$chi"
#   ../../input.py --chips 1.50 --chibs "$chi" --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --polar --bmon --phi 0.4
#   ../../stitch.py --atoms ../traj.lammpstrj
#   #~/Downloads/cuda-tild/gpu-tild
#   cd ../
#   done
# cd ../

# mkdir non-polar
# cd non-polar

#   ../input.py --chips 1.00 --chibs 2.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --midpush 0.05 --max_step 450001 --log_freq 1000 --traj_freq 10000 --phi 0.4 --bmon
#   ../stitch.py
#   # ../stitch.py --atoms traj.lammpstrj
#   ~/Downloads/cuda-tild/gpu-tild

#   for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
#   mkdir "$chi"
#   cd "$chi"
#   ../../input.py --chips 1.50 --chibs "$chi" --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --bmon --phi 0.4
#   ../../stitch.py --atoms ../traj.lammpstrj
#   #~/Downloads/cuda-tild/gpu-tild
#   cd ../
#   done
# cd ../



# rm -r AB-2d-salt-Jul-10-eq-null
# mkdir AB-2d-salt-Jul-10-eq-null
# cd AB-2d-salt-Jul-10-eq-null

# # touch info.txt
# # echo $'ABC 2d system, test the role of the parameter chi\n' >> info.txt
# # echo $'No salt added to the system\n' >> info.txt

# # #echo $'Folder polar: vary the chi_ps parameter\n' >> info.txt
# # #echo $'Folder non-polar: vary the chi_ps parameter\n' >> info.txt
# # #echo $'Folder counter-ions: chi_ps = 1.0, chi_pi varied\n' >> info.txt

# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
# cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py

# # echo $'\n\n' >> info.txt

# mkdir polar
# cd polar

#   for salt in 0.00;do #0.005 0.01 0.05 0.10 0.15 0.2 0.25 0.30; do
#   mkdir "$salt"
#   cd "$salt"
#   ../../input.py --chips 1.50 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --polar --phi 0.4 --salt "$salt" --midpush 0.05 --max_step 450001 --log_freq 1000 --traj_freq 10000
#   ../../stitch.py
#   cd ../
#   done
# cd ../

# mkdir non-polar#!/usr/bin/env python3
# cd non-polar

#   for salt in 0.00;do 
#   mkdir "$salt"
#   cd "$salt"
#   ../../input.py --chips 1.50 --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --salt "$salt" --midpush 0.05 --max_step 450001 --log_freq 1000 --traj_freq 10000
#   ../../stitch.py
#   cd ../
#   done
# cd ../



#   cd ../
#   done
# cd ../
#!/bin/bash

#for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
#for salt in 0.005 0.01 0.05 0.10 0.15 0.2 0.25 0.30; do
  # ~/Downloads/cuda-tild/gpu-tild

clear

cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/spice-up.py spice.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/write-data.bash write.bash


  mkdir ABC-100
  cd ABC-100

  ../input.py --chips 1.00 --chibs 1.50 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --midpush 0.05 --max_step 750001 --log_freq 5000 --traj_freq 50000 -N 33 --bmon --polar
  ../stitch.py
    ~/Downloads/cuda-tild/gpu-tild

  for chi in 0.00 0.25 0.50 1.00 1.25 1.50 2.00 2.50 3.00 3.50; do
  mkdir "$chi"
  cd "$chi"
  ../input.py --chips 0.25 --chibs "$chi" --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 -N 33 --bmon --polar
  ../../stitch.py --atoms ../traj.lammpstrj
  cd ../
    done
cd ../

mkdir AB-100
  cd AB-100

  ../input.py --chips 1.00 --chibs 1.50 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --midpush 0.05 --max_step 750001 --log_freq 5000 --traj_freq 50000 -N 50 --polar
  ../stitch.py
    ~/Downloads/cuda-tild/gpu-tild

  for chi in 0.00 0.10 0.25 0.50 1.00; do
  mkdir "$chi"
  cd "$chi"
  ../../input.py --chips "$chi" --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 -B -N 50 --polar
  ../../stitch.py --atoms ../traj.lammpstrj
  cd ../
    done
cd ../


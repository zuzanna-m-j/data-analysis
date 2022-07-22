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


  mkdir lamella
  cd lamella

  ../input.py --chips 1.00 --chibs 1.50 --lx 70 --ly 400 --lz 5 --Nx 151 --Ny 847 --Nz 1 --dim 2 --small --phi 0.4 --midpush 0.05 --max_step 750001 --log_freq 5000 --traj_freq 50000 -N 25 --bmon --polar
  ../stitch.py
    ~/Downloads/cuda-tild/gpu-tild

  for wlc in 0.25 0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00; do
  mkdir "$wlc"
  cd "$wlc"
  ../../input.py --chips 1.25 --chibs 0.30 --lx 70 --ly 400 --lz 5 --Nx 151 --Ny 847 --Nz 1 --dim 2 --small --phi 0.4 -N 25 --bmon --polar --wlc "$wlc"
  ../../stitch.py --atoms ../traj.lammpstrj
  cd ../
    done
cd ../

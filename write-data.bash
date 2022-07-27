#!/bin/bash

#for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
#for salt in 0.005 0.01 0.05 0.10 0.15 0.2 0.25 0.30; do
  # ~/Downloads/cuda-tild/gpu-tild


cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/spice-up.py spice.py


  cd non-polar
  for chi in 0.00 0.05 0.10 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
mkdir "$chi"
cd "$chi"
 ../../input.py --chips "$chi" --chibs 0.0 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --phi 0.4 -N 25
 ../../stitch.py --atoms ../positions.lammpstrj
cd ../
done

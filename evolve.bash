!/bin/bash

clear
echo "Evolving files..."


cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/gen-input.py input.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py stitch.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-eq.py stitch-eq.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/spice-up.py spice.py
cp /home/jello/Downloads/cuda-bd-master/research/data-analysis/write-data.bash write.bash


# cd polar
# for i in */; do
#  cd $i
#  chi=$(echo $i | sed "s/.$//")
#   rm *.tec
#   rm *.bin
#   rm positions*
#   python3 ../../input.py --chips "$chi" --chibs 0.00 --lx 50 --ly 200 --lz 5 --Nx 111 --Ny 451 --Nz 1 --dim 2 --small --phi 0.4 --polar --max_step 2000001 --traj_freq 20000
#   python3 ../../stitch.py --atoms traj.lammpstrj
#     cd ../
#     done
# cd ../


cd polar
  for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
  mkdir "$chi"
  cd "$chi"
# for i in */; do
#  cd $i
#  chi=$(echo $i | sed "s/.$//")
#  echo "$chi"
  # rm *.tec
  # rm *.bin
  # rm positions*
  
  python3 ../../input.py --chips "$chi" --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --phi 0.3 --max_step 4000001 --polar
  python3 ../../stitch.py --atoms ../positions.lammpstrj
    cd ../
    done
cd ../

cd non-polar
  for chi in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00 3.00; do
  mkdir "$chi"
  cd "$chi"
# for i in */; do
#  cd $i
#  chi=$(echo $i | sed "s/.$//")
#  echo "$chi"
  # rm *.tec
  # rm *.bin
  # rm positions*

  python3 ../../input.py --chips "$chi" --chibs 0.00 --lx 20 --ly 20 --lz 125 --Nx 51 --Ny 51 --Nz 319 --dim 3 --small --phi 0.3 --max_step 4000001
  python3 ../../stitch.py --atoms ../positions.lammpstrj
    cd ../
    done
cd ../


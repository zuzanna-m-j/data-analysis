#!/bin/bash

for i in */; do
 cd $i
shopt -s nullglob
 numdirs=(*/)
 numdirs=${#numdirs[@]}
 if [ $numdirs == 0 ]
 then
   rm input.data
   rm data.dat
   rm positions.bin
   rm grid_densities.bin
   rm traj.lammpstrj
   /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py --atoms positions.lammpstrj
   rm positions.lammpstrj
cd ../

 else
 for k in */
 do
  cd $k
  numfiles=(*)
  numfiles=${#numfiles[@]}
  if [ $numfiles > 0 ]
  then
   rm input.data
   rm data.dat
   rm positions.bin
   rm grid_densities.bin
   rm traj.lammpstrj
   /home/jello/Downloads/cuda-bd-master/research/data-analysis/stitch-input.py --atoms positions.lammpstrj
   rm positions.lammpstrj
  fi
  cd ../
 done
 cd ../
fi
done

# #!/bin/sh
# jobs=0

# for i in */; do
#  cd $i
# shopt -s nullglob
#  numdirs=(*/)
#  numdirs=${#numdirs[@]}
#  if [ $numdirs == 0 ]
#  then
#   numfiles=(*)
#   numfiles=${#numfiles[@]}
#   if [ $numfiles > 0 ]
#   then
#    rm input.data
#    rm data.dat
#    rm positions.bin
#    rm grid_densities.bin
#    rm traj.lammpstrj
#    stitch.py --atoms positions.lammpstrj
#    rm positions.lammpstrj
#   fi
#  cd ../


#  else
#  for k in */
#  do
#   cd $k
#   numfiles=(*)
#   numfiles=${#numfiles[@]}
#   if [ $numfiles > 0 ]
#   then
#    rm input.data
#    rm data.dat
#    rm positions.bin
#    rm grid_densities.bin
#    rm traj.lammpstrj
#    stitch.py --atoms positions.lammpstrj
#    rm positions.lammpstrj
#   fi
#   cd ../
#  done
#  cd ../
# fi
# done

#    rm input.data input1.data
#    rm data.dat data1.dat
#    rm positions.bin positions1.bin
#    rm grid_densities.bin grid_densities1.bin
#    rm traj.lammpstrj traj1.lammpstrj
#    stitch.py --atoms positions.lammpstrj
#    rm positions.lammpstrj



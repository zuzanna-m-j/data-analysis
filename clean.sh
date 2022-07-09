#!/bin/sh

for i in */; do
 cd $i
shopt -s nullglob
 numdirs=(*/)
 numdirs=${#numdirs[@]}
 if [ $numdirs == 0 ]
 then
  numfiles=(*)
  numfiles=${#numfiles[@]}
  if [ $numfiles > 0 ]
  then
  rm data1.dat grid_densities1.bin positions1.bin traj1.lammpstrj
  fi
 cd ../

 else
 for k in */
 do
  cd $k
  numfiles=(*)
  numfiles=${#numfiles[@]}
  if [ $numfiles > 0 ]
  then
  rm data1.dat grid_densities1.bin positions1.bin traj1.lammpstrj
  fi
  cd ../
 done
 cd ../
fi
done
echo "Clean up completed"

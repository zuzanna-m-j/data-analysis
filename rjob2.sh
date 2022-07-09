#!/bin/sh
jobs=0

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

   mv head.data head1.data
   mv input.data input1.data
   mv input input1
   mv input2 input
   mv data.dat data1.dat
   mv positions.bin positions1.bin
   mv grid_densities.bin grid_densities1.bin
   mv traj.lammpstrj traj1.lammpstrj
   stitch.py --atoms traj1.lammpstrj
   job=$(echo $i | sed "s/.$//")
   sbatch -J $job submit1-gpu;
   sleep 0.05;
   let jobs++;
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
   mv input.data input1.data
   mv data.dat data1.dat
   mv positions.bin positions1.bin
   mv grid_densities.bin grid_densities1.bin
   mv traj.lammpstrj traj1.lammpstrj
   stitch.py --atoms traj1.lammpstrj
  fi
  cd ../
 done
 cd ../
fi
done

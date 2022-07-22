#!/bin/bash

for i in */; do
 cd $i
shopt -s nullglob
 numdirs=(*/)
 numdirs=${#numdirs[@]}
 if [ $numdirs == 0 ]
 then

/home/jello/Downloads/cuda-bd-master/research/data-analysis/get-lines.py
/home/jello/results/dump-grid-densities grid_densities.bin density  $(cat skip_that)
cd ../

 else
 for k in */
 do
  cd $k
  numfiles=(*)
  numfiles=${#numfiles[@]}
  if [ $numfiles > 0 ]
  then

/home/jello/Downloads/cuda-bd-master/research/data-analysis/get-lines.py
/home/jello/results/dump-grid-densities grid_densities.bin density  $(cat skip_that)

  fi
  cd ../
 done
 cd ../
fi
done


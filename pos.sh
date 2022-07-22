#!/bin/bash

while getopts 'vf' OPTION; do
  case "$OPTION" in
    v)
    for i in */; do
    cd $i
      shopt -s nullglob
      numdirs=(*/)
      numdirs=${#numdirs[@]}
      if [ $numdirs == 0 ]
      then

      # /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj
      /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
      cd ../

      else
      for k in */
      do
        cd $k
        numfiles=(*)
        numfiles=${#numfiles[@]}
        if [ $numfiles > 0 ]
        then
      # /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj
      /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
        fi
        cd ../
      done
      cd ../
      fi
      done
    ;;

    f)
    for i in */; do
    cd $i
      shopt -s nullglob
      numdirs=(*/)
      numdirs=${#numdirs[@]}
      if [ $numdirs == 0 ]
      then
      /home/jello/Downloads/cuda-bd-master/research/data-analysis/get-lines.py
      /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj $(cat skip_that)
      /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
      cd ../

      else
      for k in */
      do
        cd $k
      /home/jello/Downloads/cuda-bd-master/research/data-analysis/get-lines.py
      /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj $(cat skip_that)
      /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
        cd ../
      done
      cd ../
      fi
      done
    ;;
  esac
done
shift "$(($OPTIND -1))"
cd ../


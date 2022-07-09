#!/bin/bash

for i in */; do
 cd $i
shopt -s nullglob
 numdirs=(*/)
 numdirs=${#numdirs[@]}
 if [ $numdirs == 0 ]
 then
name=$(echo $i| sed "s/.$//")
while getopts 'av:' OPTION; do
  case "$OPTION" in
    a)
    /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj
    ;;
    v)
    /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
    ;;
    ?)
    exit 1
    ;;
  esac
done
shift "$(($OPTIND -1))"
cd ../

 else
 for k in */
 do
  cd $k
  numfiles=(*)
  numfiles=${#numfiles[@]}
  if [ $numfiles > 0 ]
  then
        while getopts 'av:' OPTION; do
        case "$OPTION" in
            a)
            /home/jello/results/dump-particle-posits positions.bin positions.lammpstrj
            ;;
            v)
            /home/jello/Downloads/ovito-basic-3.7.4-x86_64/bin/ovito positions.lammpstrj
            ;;
            ?)
            exit 1
            ;;
        esac
        done
  fi
  cd ../
 done
 cd ../
fi
done


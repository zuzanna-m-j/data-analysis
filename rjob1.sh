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
   cp /home/zuzannaj/submit1-gpu .
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
   cp /home/zuzannaj/submit1-gpu .
   job=$(echo $k| sed "s/.$//")
   sbatch -J $job submit1-gpu;
   sleep 0.05;
   let jobs++;
  fi
  cd ../
 done
 cd ../
fi
done
echo "Zuzanna, now you have $jobs running"
sleep 20;
squeue -u zuzannaj

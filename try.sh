#!/bin/bash


Total_folder=1

for i in $(seq 1 $Total_folder)
do
    #cp jobsubmission.sh ${i}/jobsubmission.sh
    #cd ${i}
    #echo "echoMe $i" >> jobsubmission_${i}.sh
    #chmod +x jobsubmission_${i}.sh
    ./jobsubmission.sh
    #echo "echoMe $i" >> jobsubmission_${i}.sh
    #lmp_serial < in.lammps
    sleep 1
    #cd ..
done

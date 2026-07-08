#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:48:35 2023

@author: rahul
"""
import numpy
import time
import random
import os
import subprocess
import shutil
import time
start_time = time.time()




start = 0
Total_file_generate = 100
MDstep = 1000
for i in range(start,Total_file_generate):
    #shutil.copy('in-fene.lammps', str(i+1)+'/in-fene.lammps')
    #change data file content
    if i == 0:
        lammps_file = open("db_in.lammps", "r")
        lammps_content = lammps_file.readlines()
        lammps_content[24] = 'reset_timestep ' + str(i*MDstep) + '\n'
        lammps_file.close()
        lammps_file_mod = open("db_in.lammps", "w")
        lammps_file_mod.writelines(lammps_content)
        lammps_file_mod.close()
        subprocess.call("./try.sh") 
        print('i value is : ', i)
    if i > 0:
        lammps_file = open("db_in.lammps", "r")
        lammps_content = lammps_file.readlines()
        lammps_content[3] = 'read_data       system.data\n'
        lammps_content[24] = 'reset_timestep ' + str(i*MDstep) + '\n'
        if (i+1)%100 ==0 :
            lammps_content[54] = 'variable        refile string "restart.data"\n'
            lammps_content[55] = 'restart         1000 ${refile}\n'
        else:
            lammps_content[54] = '#variable        refile string "restart.data"\n'
            lammps_content[55] = '#restart         1000 ${refile}\n'
        if i == 4000:
            lammps_content[41] = 'dump            10 all custom 1000 traj-4000-6000.lammpstrj id mol type x y z ix iy iz\n'
            lammps_content[42] = 'dump_modify     10 sort id append yes\n'
        if i+1 == 6000:
            lammps_content[41] = 'dump            10 all custom 1000 traj-6000-8000.lammpstrj id mol type x y z ix iy iz\n'
            lammps_content[42] = 'dump_modify     10 sort id append yes\n'
        if i+1 == 8000:
            lammps_content[41] = 'dump            10 all custom 1000 traj-8000-10000.lammpstrj id mol type x y z ix iy iz\n'
            lammps_content[42] = 'dump_modify     10 sort id append yes\n'
            
        lammps_file.close()
        lammps_file_mod = open("db_in.lammps", "w")
        lammps_file_mod.writelines(lammps_content)
        lammps_file_mod.close()
        system_file = open("system.data", "r")
        list_of_content = system_file.readlines()
        list_of_content[16] = '\n'
        list_of_content[18] = '\n'
        list_of_content[19] = '\n'
        list_of_content[21] = '\n'
        list_of_content[23] = '\n'
        list_of_content[24] = '\n'
        system_file.close()
        system_file_mod = open("system.data", "w")
        system_file_mod.writelines(list_of_content)
        system_file_mod.close()
        subprocess.call("./try.sh")
        print('i value is : ', i)
    # else:  subprocess.call("./try.sh") 
    subprocess.run(["python", "MC.py"])
    print('i value is : ', i)
end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.6f} seconds")   

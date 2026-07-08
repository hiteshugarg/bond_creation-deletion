#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 16:27:46 2023

@author: rahul
"""

import math
import random
import time
import numpy as np

OUTPUT_FILE = "MCoutput.txt"

start_time = time.time()

def _pbc_distance(x1, y1, z1, x2, y2, z2, box_dims):
    xdis = x1 - x2
    xdis = xdis - box_dims[0] * round(xdis / box_dims[0])
    ydis = y1 - y2
    ydis = ydis - box_dims[1] * round(ydis / box_dims[1])
    zdis = z1 - z2
    zdis = zdis - box_dims[2] * round(zdis / box_dims[2])
    return math.sqrt(xdis * xdis + ydis * ydis + zdis * zdis)


def _log(line):
    with open(OUTPUT_FILE, "a+", encoding="utf-8") as file:
        file.write(line)


def _choose_weighted(candidates, weights):
    if not candidates:
        return None
    assert len(candidates) == len(weights)
    total = float(sum(weights))
    if total <= 0.0:
        return None

    pick = random.random() * total
    cumulative = 0.0
    for candidate, weight in zip(candidates, weights):
        cumulative += float(weight)
        if pick <= cumulative:
            return int(candidate)
    return int(candidates[-1])


# Create Neighbourlist and distancelist in python
def build_verlet_list2(x,y,z, box_size, cutoff,chain_length): #verlet neighbour list with 1st method
    verlet_list = [[] for _ in range(N_tot)]
    distance_list = [[] for _ in range(N_tot)]
    for i in range(N_tot-1):
        if (i+1) % chain_length ==0:
            for j in range(i + 1, N_tot):
                xdis = x[i] - x[j]
                xdis = xdis - box_size[0][0]*round(xdis/box_size[0][0])
                ydis = y[i] - y[j]
                ydis = ydis - box_size[0][1]*round(ydis/box_size[0][1])
                zdis = z[i] - z[j]
                zdis = zdis - box_size[0][2]*round(zdis/box_size[0][2])
                distance = math.sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
                if distance < cutoff:
                    #verlet_list[atom_id[i]-1].append(atom_id[j]-1)
                    #verlet_list[atom_id[j]-1].append(atom_id[i]-1)
                    verlet_list[i].append(j)
                    verlet_list[j].append(i)
                    distance_list[i].append(distance)
                    distance_list[j].append(distance)
        else:
            for j in range(i + 2, N_tot):
                xdis = x[i] - x[j]
                xdis = xdis - box_size[0][0]*round(xdis/box_size[0][0])
                ydis = y[i] - y[j]
                ydis = ydis - box_size[0][1]*round(ydis/box_size[0][1])
                zdis = z[i] - z[j]
                zdis = zdis - box_size[0][2]*round(zdis/box_size[0][2])
                distance = math.sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
                if distance < cutoff:
                    #verlet_list[atom_id[i]-1].append(atom_id[j]-1)
                    #verlet_list[atom_id[j]-1].append(atom_id[i]-1)
                    verlet_list[i].append(j)
                    verlet_list[j].append(i)
                    distance_list[i].append(distance)
                    distance_list[j].append(distance)
                
    return verlet_list , distance_list

def build_verlet_list(x,y,z, box_size, cutoff,chain_length): #verlet neighbour list with 2nd method
    verlet_list = [[] for _ in range(N_tot)]
    distance_list = [[] for _ in range(N_tot)]
    for i in range(N_tot-1):
        if atom_id[i] % chain_length == 0:
            for j in range(i + 1, N_tot):
                if atom_id[j] != atom_id[i]-1:
                    xdis = x[i] - x[j]
                    xdis = xdis - box_size[0][0]*round(xdis/box_size[0][0])
                    ydis = y[i] - y[j]
                    ydis = ydis - box_size[0][1]*round(ydis/box_size[0][1])
                    zdis = z[i] - z[j]
                    zdis = zdis - box_size[0][2]*round(zdis/box_size[0][2])
                    distance = math.sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
                    if distance < cutoff:
                        #verlet_list[atom_id[i]-1].append(atom_id[j]-1)
                        #verlet_list[atom_id[j]-1].append(atom_id[i]-1)
                        verlet_list[i].append(j)
                        verlet_list[j].append(i)
                        distance_list[i].append(distance)
                        distance_list[j].append(distance)
                #else: print(atom_id[j], atom_id[i])        
        else:
            for j in range(i + 1, N_tot):
                if atom_id[j] != atom_id[i]-1 and atom_id[j] != atom_id[i]+1:
                    xdis = x[i] - x[j]
                    xdis = xdis - box_size[0][0]*round(xdis/box_size[0][0])
                    ydis = y[i] - y[j]
                    ydis = ydis - box_size[0][1]*round(ydis/box_size[0][1])
                    zdis = z[i] - z[j]
                    zdis = zdis - box_size[0][2]*round(zdis/box_size[0][2])
                    distance = math.sqrt(xdis*xdis+ydis*ydis+zdis*zdis)
                    if distance < cutoff:
                        #verlet_list[atom_id[i]-1].append(atom_id[j]-1)
                        #verlet_list[atom_id[j]-1].append(atom_id[i]-1)
                        verlet_list[i].append(j)
                        verlet_list[j].append(i)
                        distance_list[i].append(distance)
                        distance_list[j].append(distance) 
                #else: print(atom_id[j], atom_id[i])        
    return verlet_list , distance_list

# Read neighbourlist and distancelist from Lammps dump file
def read_pair_distance(filename,verlet_list_test,distance_list_test):
    i = 0
    fpair = open(filename , 'r')
    while True:
        s = fpair.readline().split()
        if not s : break
        if i == 3:
            entry = int(s[0])
            break
        i = i + 1
    #print(entry)  
    fpair.close()
    fpair = open(filename , 'r') 
    i = 0
    while True:
        s = fpair.readline().split()
        if not s : break
        if i%(entry+9) == 0 :
            j = 0
        if j==1:
            timestep = int(s[0])
        if j > 8:
            atom1 = int(s[1])
            atom2 = int(s[2])
            dis12 = float(s[3])
            if dis12 < cutoff:
                if bead_id[atom1-1]*bead_id[atom2-1] == 2:
                        verlet_list_test[atom1-1].append(atom2-1)
                        verlet_list_test[atom2-1].append(atom1-1)
                        distance_list_test[atom1-1].append(dis12)
                        distance_list_test[atom2-1].append(dis12)
        j = j + 1
        i = i + 1
    fpair.close()
    return timestep,entry,verlet_list_test,distance_list_test


def read_lammps_traj(filename,natoms):
    atom_id.clear()
    mol_id.clear()
    bead_id.clear()
    x.clear()
    y.clear()
    z.clear()

    i = 0
    with open(filename, "r", encoding="utf-8") as file:
        while True:
            s = file.readline().split()
            if not s:
                break
            if i % (natoms + 9) == 0:
                j = 0
            if j == 1:
                timestep = int(s[0])
            if j > 8:
                atom_id.append(int(s[0]))
                mol_id.append(int(s[1]))
                bead_id.append(int(s[2]))
                x.append(float(s[3]))
                y.append(float(s[4]))
                z.append(float(s[5]))
            j = j + 1
            i = i + 1

    return timestep

def read_lammps_data_file(filename):
    box_dims_list = []
    # open the file for reading
    with open(filename, 'r') as file:
        lines = file.readlines()
        natoms = int(lines[2].split()[0])
        num_crosslink = np.zeros(natoms)
        bonds  = int(lines[4].split()[0])
        Lx = float(lines[7].split()[1])-float(lines[7].split()[0])
        Ly = float(lines[8].split()[1])-float(lines[8].split()[0])
        Lz = float(lines[9].split()[1])-float(lines[9].split()[0])
        box_dims_list.append([Lx,Ly,Lz])
        #found_atoms = False
        #found_bonds = False
        #atom_lines = []
        #bond_lines = []
        #timestep = None
        crosslink_monomer = []
        
        # iterate over each line in the file
        i = 0
        while i < len(lines):
            timestep = lines[0].split()[11]
            #print(timestep)
            # if "Atoms # molecular" in lines[i]:
            #     #atom_id = []
            #     #mol_id = []
            #     #bead_id = []
            #     #x = []
            #     #y = []
            #     #z = []
            #     for j in range(natoms):
            #         atom_id.append(int(lines[i+j+2].split()[0]))
            #         mol_id.append(int(lines[i+j+2].split()[1]))
            #         bead_id.append(int(lines[i+j+2].split()[2]))
            #         x.append(float(lines[i+j+2].split()[3]))
            #         #x[atom_id[j]-1] = float(lines[i+j+2].split()[3])
            #         y.append(float(lines[i+j+2].split()[4]))
            #         #y[atom_id[j]-1] = float(lines[i+j+2].split()[4])
            #         z.append(float(lines[i+j+2].split()[5]))
            #         #z[atom_id[j]-1] = float(lines[i+j+2].split()[5])
            #         #print(atom_id[j],bead_id[j],x[j],y[j],z[j])
                    
            if "Bonds" in lines[i]:  
                #bond_array = np.zeros((bonds,4))
                #bond_id = []
                #bond_1  = []
                #bond_2 = []
                for j in range(bonds):
                    bond_id = int(lines[i+j+2].split()[1])
                    bond_1 = int(lines[i+j+2].split()[2])
                    bond_2 = int(lines[i+j+2].split()[3])
                    if bond_id ==1:
                        Nbond_id.append(bond_id)
                        Nbond_1.append(bond_1)
                        Nbond_2.append(bond_2)
                        num_crosslink[bond_1-1] += 0
                        num_crosslink[bond_2-1] += 0
                    if bond_id ==2 or bond_id ==3:    
                        #if bond_1 not in crosslink_monomer: 
                        crosslink_monomer.append(bond_1)
                        num_crosslink[bond_1-1] += 1
                        #print(num_crosslink[bond_1-1])
                        #if bond_2 not in crosslink_monomer: 
                        crosslink_monomer.append(bond_2)
                        num_crosslink[bond_2-1] += 1
                        crosslink_bond_id.append(bond_id)
                        crosslink_bond1.append(bond_1)    
                        crosslink_bond2.append(bond_2)      
            i += 1
        
    
    
    
    return timestep, natoms,bonds, box_dims_list, crosslink_monomer,num_crosslink 


crosslink_bead_type = [2]
atom_id = []
mol_id = []
bead_id = []
x = []
y = []
z = []

Nbond_id = []
Nbond_1  = []
Nbond_2 = []
crosslink_bond_id = []
crosslink_bond1 = []
crosslink_bond2 = []
filea = "system.data" #Read Bond Information

#time_aa_ab_bb, count2_aa_ab_bb, count3_aa_ab_bb = read_lammps_data_file(file)
time_aa_ab_bb, N_tot,Nbonds,box_dims_list, crosslink_monomer,Num_crosslink = read_lammps_data_file(filea)

#atom_id.sort()  # Sort the list in descending order
U_bond = len(Nbond_id)

M = N_tot - len(crosslink_monomer)

Nb = len(crosslink_monomer)/2


fileb ="traj-current.lammpstrj" #Read Coordinate Information


time_aa_ab_bb = read_lammps_traj(fileb,N_tot)

Reactive_beads =  [atom_id[i] for i in range(len(bead_id)) if bead_id[i] in crosslink_bead_type]
Polymers = N_tot-len(Reactive_beads)
chain_length = int(Polymers/(Polymers-U_bond))
print('Current_timestep: ' , time_aa_ab_bb)
print('Total Number of atoms:', N_tot)
print('Total Number of Polymers:', Polymers)
print('Total normal bonds: ' ,U_bond)
print('Chain length in no crosslink state is:', chain_length)
print('Number of crosslink bonds:', Nb)
cutoff = 1.49
   
# These are written for crosscheck with neighbourlist value in python , C++ and Lammps
# verlet_list1, distance_list1 = build_verlet_list2(x,y,z,bead_id, box_dims_list, cutoff,chain_length)

#%%%
# Create Neighbourlist and distancelist in C++ and read them
# Compile the C++ program (if needed)
# subprocess.run(['g++', 'neighbourlist.cpp', '-o', 'neighbourlist'])




# result = subprocess.run(
#     ['./neighbourlist', str(N_tot), str(Polymers),str(U_bond),str(cutoff),str(crosslink_bead_type)],
#     capture_output=True,
#     text=True
# )

# print(result.stdout.strip()) 


# file_neighbourlist = open('neighbourlist.txt' , 'r')
# verlet_list2 = [[] for _ in range(N_tot)]
# i = 0
# while True:
#     line = file_neighbourlist.readline()
#     s = line.split()
#     # if not s.split: 
#     for j in range(len(s)):
#         verlet_list2[i].append(int(s[j]))
#     if not line: break    
#     i += 1    

# file_distancelist = open('distancelist.txt' , 'r')
# distance_list2 = [[] for _ in range(N_tot)]
# i = 0
# while True:
#     line = file_distancelist.readline()
#     s = line.split()
#     # if not s.split: 
#     for j in range(len(s)):
#         distance_list2[i].append(float(s[j]))
#     if not line: break    
#     i += 1    
#%%%

verlet_list = [[] for _ in range(N_tot)]
distance_list = [[] for _ in range(N_tot)]
filepairdis = "pairdistance.dump"            

time_pair_file,entry,verlet_list, distance_list= read_pair_distance(filepairdis,verlet_list,distance_list)

print('Current_timestep in pair file: ' , time_pair_file)
print('Total entry in pair file: ' , entry)

# Incorporation of crosslinked monomer in the neighbourlist as they were not present in the list of pair distance
for i in range(len(crosslink_bond1)):
    pair1 = crosslink_bond1[i]
    pair2 = crosslink_bond2[i]
    verlet_list[pair1-1].append(pair2-1)
    verlet_list[pair2-1].append(pair1-1)
    bondeddis = _pbc_distance(
        x[pair1 - 1],
        y[pair1 - 1],
        z[pair1 - 1],
        x[pair2 - 1],
        y[pair2 - 1],
        z[pair2 - 1],
        box_dims_list[0],
    )
    assert(bondeddis < cutoff)
    distance_list[pair1-1].append(bondeddis)
    distance_list[pair2-1].append(bondeddis)

# verlet_list3   = [sorted(row) for row in verlet_list] 
#%%
# if verlet_list2 == verlet_list3 :
#     print("Both lists are equal.")
# else:
#     print("Lists are not equal.")


# if distance_list2 == distance_list :
#     print("Both lists are equal.")
# else:
#     print("Lists are not equal.")
#%%
_log("##----------------------------##\n")
_log("timestep = %s\n" % (str(time_aa_ab_bb)))
#start MC
MCstep = 200
Pf = Pr = 0.5
mu = 80 # This is the dynamic bond chemical potential, It will dictate the crosslink percentage

T = 1
K = 30.0
sigma = 1.0
Rc = 1.5*sigma
eps = 1.0
print('Chemical Potential:', mu)
print('Temperature:', T)
print('paair interaction strength:', eps)
print('particle diameter:', sigma)
print('Fene Bond Strength:', K)
print('Fene Bond max Length:', Rc)
def FENE(r,K,Rc,eps,T):
    U = -0.5*K*Rc*Rc*math.log(1-(r/Rc)*(r/Rc)) + 4*eps*((sigma/r)**12-(sigma/r)**6) + eps
    #print(U)
    U = math.exp(-U/T)
    return U

def chooseAtom(iatom,Vlist,Blist): #choose jth atom from ith neighbour list
    return _choose_weighted(Vlist, Blist)

#start MC for bond creation
def Bond_creation(Nb,M):
    #while True:
    Ith = random.choice(Reactive_beads) # The probability to choose ith particle randomly
    Ithindex = atom_id.index(Ith) 
    #if Ith not in crosslink_monomer: 
        #break
    
    sample = []  
    boltzman_factor_list = []
    if Num_crosslink[Ithindex] < 2:
        if verlet_list[Ithindex]:
            Wi = 0 # initialize the Wi 
            for neighbor_index, neighbor_distance in zip(
                verlet_list[Ithindex], distance_list[Ithindex]
            ):
                # eligible partner must have < 1 existing crosslink
                if Num_crosslink[neighbor_index] < 1:
                    wth = FENE(neighbor_distance, K, Rc, eps, T)
                    Wi = Wi + wth
                    boltzman_factor_list.append(wth)
                    sample.append(neighbor_index)  # possible j who can form new crosslink
        else:
            _log("%s no neighbors\n" % (str(Ith)))
    else: 
        _log("%s max crosslinks reached\n" % (str(Ith)))
    if sample:
        Jthindex = chooseAtom(Ithindex,sample,boltzman_factor_list) 
        if Jthindex is None:
            return Nb, M
        Jth = atom_id[Jthindex]
        #print(Jthindex,verlet_list[Ithindex])
        
        Wj = 0  # initialize the Wj
        for j in range(len(verlet_list[Jthindex])):
            #if atom_id[verlet_list[Jthindex][j]] not in crosslink_monomer:
            if Num_crosslink[verlet_list[Jthindex][j]] < 2:
                Wj = Wj+ FENE(distance_list[Jthindex][j],K,Rc,eps,T)  #summing over all possible bonding partner  
                
            
        Wij = (Wi*Wj)/(Wi+Wj) 
        Z = math.exp(mu/T)
        Bc =  (Pr/Pf)*(M/(Nb+1))*Z*Wij 

        randc = random.random() #random number for creating bond
        print(randc, Bc)
        if randc < Bc :
            if bead_id[Ithindex] == bead_id[Jthindex]:
                crosslink_bond_id.append(3)
                crosslink_bond1.append(atom_id[Ithindex])
                crosslink_bond2.append(atom_id[Jthindex])
            else:
                crosslink_bond_id.append(2)
                crosslink_bond1.append(atom_id[Ithindex])
                crosslink_bond2.append(atom_id[Jthindex])
            #if atom_id[Ithindex] not in crosslink_monomer:     
            crosslink_monomer.append(atom_id[Ithindex])
            #if atom_id[Jthindex] not in crosslink_monomer:     
            crosslink_monomer.append(atom_id[Jthindex])
            count_of_numberI = crosslink_monomer.count(atom_id[Ithindex])
            count_of_numberJ = crosslink_monomer.count(atom_id[Jthindex])
            #print(count_of_numberI,count_of_numberJ)
            assert(count_of_numberI <= 2)
            assert(count_of_numberJ < 2)
            assert(Num_crosslink[Ithindex] >= 0 and Num_crosslink[Ithindex] < 2)
            assert(Num_crosslink[Jthindex] >= 0 and Num_crosslink[Jthindex] < 1)
            Num_crosslink[Ithindex] += 1
            Num_crosslink[Jthindex] += 1 
            Nb = Nb + 1    
            if Num_crosslink[Ithindex] == 2:
                M = M - 2
            else:
                M = M - 1
            _log("create bond %s %s\n" % (str(Ith), str(Jth)))
        #break
        else: 
            _log("not created bond %s %s\n" % (str(Ith), str(Jth)))
    return Nb, M   


#start MC for bond breaking        
def Bond_breaking(Nb,M):
    if crosslink_bond1:
        Ith = random.choice(crosslink_bond1) #choose a bond i-j ith one
        
        Ithindex = crosslink_bond1.index(Ith) #choose a bond i-j jth one
        
        Jthindex = Ithindex
        Jth = crosslink_bond2[Jthindex]
        
        atom_id_pos1 = atom_id.index(Ith) 
        atom_id_pos2 = atom_id.index(Jth) 
        
        Wi = 0 # initialize the Wi
        
        #for i in range(len(verlet_list[atom_id_pos1])): #summing over all possible bonding partner
             #if atom_id[verlet_list[atom_id_pos1][i]] not in crosslink_monomer :
                 #Wi = Wi+ FENE(distance_list[atom_id_pos1][i],K,Rc,eps,T)
        
        Wj = 0 # initialize the Wj
        #for j in range(len(verlet_list[atom_id_pos2])): #summing over all possible bonding partner
             #if atom_id[verlet_list[atom_id_pos2][j]] not in crosslink_monomer :       
                 #Wj = Wj+ FENE(distance_list[atom_id_pos2][j],K,Rc,eps,T)
        
        # Bonded atoms contribution    
        bondeddis = _pbc_distance(
            x[atom_id_pos1],
            y[atom_id_pos1],
            z[atom_id_pos1],
            x[atom_id_pos2],
            y[atom_id_pos2],
            z[atom_id_pos2],
            box_dims_list[0],
        )
        assert(bondeddis < cutoff)
        bondedwij = FENE(bondeddis,K,Rc,eps,T)
        Wi = Wi + bondedwij
        Wj = Wj + bondedwij
        
        Wij = (Wi*Wj)/(Wi+Wj) 
        Z = math.exp(mu/T)
        Bb =  (Pf/Pr)*((Nb)/M/Z/Wij) ### Note: we use Nb instead of Nb + 1 as we are breaking bond.
        
        randb = random.random() #random number for breaking bond
        if randb < Bb :
            del crosslink_bond_id[Ithindex]
            del crosslink_bond1[Ithindex]
            del crosslink_bond2[Jthindex]
            crosslink_monomer.remove(Ith)
            crosslink_monomer.remove(Jth)  
            #print('Number of crosslink',Num_crosslink[atom_id_pos1])
            assert(Num_crosslink[atom_id_pos1] >= 1 and Num_crosslink[atom_id_pos1] <= 2)
            assert(Num_crosslink[atom_id_pos2] == 1)
            Num_crosslink[atom_id_pos1] -= 1
            Num_crosslink[atom_id_pos2] -= 1 
            Nb = Nb - 1
            M = M + 2
            _log("delete bond %s %s\n" % (str(Ith), str(Jth)))
        else: 
            _log("not delete bond %s %s\n" % (str(Ith), str(Jth)))
    return Nb, M


breaking_step = 0
creating_step = 0
for MC in range(MCstep):
    #random_break_create = random.random()
    #if random_break_create > Pf:
        Nb,M = Bond_creation(Nb,M)
        creating_step += 1 
    #else:
        #Nb,M = Bond_breaking(Nb,M)   
        #breaking_step += 1
    
with open(filea, "r", encoding="utf-8") as system_file_bond:
    bond_contents = system_file_bond.readlines()
del bond_contents[-Nbonds:]
bond_contents[4] = str(int(Nb+U_bond))+ '  bonds\n' 
for number in range(U_bond):
    bond_contents.append(str(number+1) + ' ' + str(Nbond_id[number]) + ' ' + str(Nbond_1[number]) + ' ' + str(Nbond_2[number]) +'\n')
for number in range(int(Nb)):
    bond_contents.append(str(number+1+U_bond)+ ' ' + str(crosslink_bond_id[number]) + ' ' + str(crosslink_bond1[number]) + ' ' + str(crosslink_bond2[number]) +'\n')
with open(filea, "w", encoding="utf-8") as system_file_bond_update:
    system_file_bond_update.writelines(bond_contents)
end_time = time.time()
elapsed_time = end_time - start_time
print('Creation attempt', creating_step/MCstep)
print('Breaking attempt', breaking_step/MCstep) 
print(f"Elapsed time: {elapsed_time:.6f} seconds")   

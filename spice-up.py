#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--head', default = "head.data")
parser.add_argument('--atoms', default = "../traj.lammpstrj")
parser.add_argument('--tail', default = "../tail.data")
parser.add_argument('--file', default = "input.data")


args = parser.parse_args()

header_file = args.head
atom_file = args.atoms
tail_file = args.tail
file_name = args.file
spice_info = np.loadtxt('spice.info', skiprows = 1)
n_salt = int(spice_info[1])
n_cat = n_salt//2
n_anion = n_salt//2
n_salt = n_cat + n_anion
n_sol = int(spice_info[0])

print(n_sol)
print(n_salt)

with open("check", 'w') as f:
    f.writelines("Hi")

with open(file=header_file,mode='r') as f:
    head = f.readlines()

with open(file=atom_file,mode='r') as f:
    atom_lines = f.readlines()

n_atoms = atom_lines[3]    

with open(file=tail_file,mode='r') as f:
    tail = f.readlines()

with open(file=file_name,mode='w') as fout:
    
    for line in head:
        fout.writelines(line)

    fout.writelines('\n')
    fout.writelines('Atoms\n')
    fout.writelines('\n')

    for atom in atom_lines[-int(n_atoms):-n_sol]:
        line = atom.split()
        if len(line) == 7:
            fout.writelines(f"{line[0]} {line[2]} {line[1]}  {line[6]}  {line[3]}  {line[4]}  {line[5]}\n")

    solvent = atom_lines[-n_sol:]

    water = solvent[:n_sol-n_salt*2]

    last_water = water[-2]
    last_water = last_water.split()

    last_mol = int(last_water[2])
    last_atom = int(last_water[0])

    atom_count = last_atom + 2

    salt = solvent[n_sol-n_salt*2::2]


    for atom in water:
        line = atom.split()
        if len(line) == 7:
            fout.writelines(f"{line[0]} {line[2]} {line[1]}  {line[6]}  {line[3]}  {line[4]}  {line[5]}\n")
    print(f"{line[0]} {line[2]} {line[1]}  {line[6]}  {line[3]}  {line[4]}  {line[5]}\n")

    for atom in salt:
        line = atom.split()
        if len(line) == 7:
            print(f"{atom_count} {line[2]} {1.000}  {line[6]}  {line[3]}  {line[4]}  {line[5]}\n")
            fout.writelines(f"{atom_count} {line[2]} {1.000}  {line[6]}  {line[3]}  {line[4]}  {line[5]}\n")
            atom_count += 1

    # for atom in salt[n_salt//2:]:
    #     line = atom.split()
    #     if len(line) == 7:
    #         fout.writelines(f"{line[0]} {6} {line[1]}  {line[6]}  {line[3]}  {line[4]}  {-1.000}\n")
    #         atom_count += 1
    #         mol_count += 1

    for line in tail:
        fout.writelines(line)
        tail_line = line.split()
        if len(tail_line) > 2:
            if int(tail_line[2]) == last_atom:
                break


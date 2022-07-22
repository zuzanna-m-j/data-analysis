#!/usr/bin/env python3

from math import floor


with open('data.dat','r') as f:
    lines = f.readlines()
    lines = lines[-1]
    line = lines.strip().split()
    steps = int(line[0])

with open("input",'r') as f:
    lines = f.readlines()
    lines = lines[4]
    line = lines.strip().split()
    traj_freq = int(line[1])

num_traj = floor(steps/traj_freq)
skipped = num_traj - 10

print(f"Total frames: {num_traj}")
print(f"Skipping: {skipped}")

with open("skip_that",'w') as f:
    f.writelines(f"{skipped}")
    
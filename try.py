#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import random

parser=argparse.ArgumentParser()
parser.add_argument('--myflag', action='store_true')
args=parser.parse_args()
myflag = args.myflag
print(myflag)

A =   [1,+1, 0]
B =   [2, 0, 0]
C =   [3,-1, 0]
W =   [4, 0, 0] 
CAT = [5,+1, 0]
ANI = [6,-1, 0]
CI =  [7, 1, 0]
D =   [8, 1, 0]

def SCD(sq):
    print(sq)
    sq = list(sq)
    print(sq)
    print(len(sq))
    s = []
    for m in sq:
        if m == "A":
            s.append(A[1])
        elif m == "C":
            s.append(C[1])
        else:
            s.append(0)
    scd = 0
    print(s)
    for m in range(1,len(s)):
        for n in range(0,m):
            print(m,n)
            scd += s[m] * s[n] * np.sqrt(m - n)
    scd = scd/len(s)
    return scd

BMON = True

N_a = 25
if BMON == True:
    N_b = 25
else:
    N_b = 0
N_c = 25
N = N_a + N_b + N_c
seq = N_a * "A" + N_b * "B" + N_c * "C"

seq_l = list(seq)
random.shuffle(seq_l)
seq = ''.join(seq_l)

SCD(seq)
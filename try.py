#!/usr/bin/env python3
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--myflag', action='store_true')
args=parser.parse_args()
myflag = args.myflag
print(myflag)
#!/usr/bin/python3

import sys
import gzip
import argparse
import numpy as np
import pandas as pd
from itertools import product


# Setup for the k-space grid
# Set up the grid for tract lengths
kgrid = np.arange(1,1000,1)
kgrid = np.append(kgrid,np.arange(1000,10000,10))
kgrid = np.append(kgrid,np.arange(10000,100001,100))

def help_fun2(start,maximum):
    return lambda x: np.minimum(np.maximum(0,x-start),maximum)

def help_fun(offset,a,b):
    return lambda x: help_fun2(offset,max(a,b))(x)-help_fun2(offset+min(a,b),max(a,b))(x)

def computeTractFunc(length,a,b):
    sa=[sum(a[:i]) for i in range(len(a))]
    sb=[sum(b[:i]) for i in range(len(b))]
    return lambda x: np.sum([help_fun(length+sa[i]+sb[j],a[i],b[j])(x)*(1-penetrance)**(i+j) for (i,j) in product(range(len(a)),range(len(b))) if length+sa[i]+sb[j]<=kgrid[-1]], axis=0 )

parser = argparse.ArgumentParser(description='calculate tract functions from tracts and writes them to the file [tracts].p[PENETRANCE].twdet.gz - The output file is suitable as input for EMmodel.py')
parser.add_argument('tracts', type=str, help='file containing the tracts')
parser.add_argument('--penetrance', type=float,default=0.5, help='probability of NCO marker becoming gene converted')
args = parser.parse_args()

tracts=pd.read_csv(args.tracts,delimiter='\t')
penetrance=args.penetrance

twdet=[]
count=0

def dist2diff(gaps):
    return [gaps[0]]+[gaps[i+1]-gaps[i]  for i in range(len(gaps)-1)]

for t in tracts.itertuples():
    lgap=[int(i) for i in t.lgap.split(',')]
    rgap=[int(i) for i in t.rgap.split(',')]
    # The gaps are distance from last gc marker to flinking markers
    # Change it to distance to adjacent markers
    lgap=dist2diff(lgap)+[kgrid[-1]+1]
    rgap=dist2diff(rgap)+[kgrid[-1]+1]

    length=t.ubd-t.lbd
    if length>=kgrid[-1]:
        print(count,":","length>",kgrid[-1],":",length)
        continue

    twdet.append(computeTractFunc(length,lgap,rgap))
    count+=1

twfile=args.tracts+".p%s.twdet.gz"%penetrance
print(count,"functions created. Output will be written to "+twfile)
twdet_grid=[tf(kgrid) for tf in twdet]
np.savetxt(twfile,twdet_grid)


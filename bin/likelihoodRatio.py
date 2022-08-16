#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
from scipy.special import xlogy
from scipy import stats
import argparse
from collections import namedtuple
from scipy.stats.distributions import chi2

def parseResult(filename):
    header=None
    like=None
    with open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if header is None: 
                header=line.strip().split()
            else:
                like=line.strip().split()
        return header,like

def getNgroups(header):
    n=0
    while True:
        if f"partition{n}" not in header:
            break
        n+=1
    return n


parser = argparse.ArgumentParser(description='Tekes multiple output files of the EM algorithm and calculates the likelihood ratio between the outcomes')
parser.add_argument('resultfile', nargs='+', help='output of the EM algorithm')
parser.add_argument('-v','--verbose',action='store_true', help='print more')
theargs=parser.parse_args()

entries=[]

for resultfile in theargs.resultfile:
    header,line=parseResult(resultfile)
    loglike=float(dict(zip(header,line))['LogLikelihood'])
    shapeparams=0
    fam=resultfile.split('/')[-1].split('.')[0]
    if fam=='geom':
        shapeparams=1
    elif fam=='nbinom':
        shapeparams=2
    ngroups=getNgroups(header)
    dof=(shapeparams+1)*ngroups-1

    if theargs.verbose:
        print(resultfile,loglike,dof,shapeparams,fam)

    entries.append([resultfile,loglike,dof])

print(*[f'{i}{j}' for i in ['1st','2nd'] for j in ['File','LogL','DegF']],"p-value",sep='\t')
for e1,e2 in zip(entries[:-1],entries[1:]):
    n1,l1,d1=e1
    n2,l2,d2=e2

    print(*e1,*e2,chi2.sf(2*(l2-l1),d2-d1),sep='\t')


    



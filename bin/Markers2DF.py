#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import argparse
from math import floor

def read_markerfile(filepath):
  '''
  Read gzipped markerfile with columns Chr and pos, tab separated.
  '''
  markers = pd.read_csv(filepath,sep='\t')
  return markers

def flatten_genome(markers,maxgap = 100000):
  '''
  Flatten the markerset into a contiguous sequence with maxgap (100kb) 
  between chromosomes and starting at zero on chromosome 1.
  The chromosomes are flattened in alphabetical order. Centromeres and 
  large caps are replaced with a gap of size maxgap (100kb).
  '''
  markers["diff"] = markers.groupby("Chr").diff()['pos'].fillna(maxgap)
  markers["diff"] = markers["diff"].apply(lambda x:min(maxgap,x))
  markers.at[0,"diff"] = 0
  markers["Pos"] = markers["diff"].cumsum().astype("int64")
  return list(markers["Pos"])


klist = [1]
klist += list(range(10,100,10))
klist += list(range(100,1000,100))
klist += list(range(1000,10000,1000))
klist += list(range(10000,100001,10000))

def dfprob(val,diffs):
  # Calculate detecion probability for the genome
  ExpVal=0
  thelist=[]

  # Slide a window of size val over region to calculate sum of expected #found gene conversion tracts
  # over all placements the gene conversion tract
  for d in diffs:
    # The window is empty skip ahead until next marker enters window 
    # this case can only occur for first marker in region
    if not thelist:
      thelist.append(val)
      continue
    minv=min(thelist)

    # Check for markers exiting window before next marker enters it
    while minv<d:
      ExpVal+=minv*(1-(1-penetrance)**len(thelist))
      thelist.remove(minv)
      if not thelist:
        break
      thelist=[v-minv for v in thelist]
      d-=minv
      minv=min(thelist)
    # Skip ahead until next marker enters window
    ExpVal+=d*(1-(1-penetrance)**len(thelist))
    thelist=[v-d for v in thelist]
    if 0 in thelist:
      thelist.remove(0)
    thelist.append(val)

  # slide window until all markers have exited the window
  while thelist:
    minv=min(thelist)
    ExpVal+=minv*(1-(1-penetrance)**len(thelist))
    thelist.remove(minv)
    thelist=[v-minv for v in thelist]

  return ExpVal


def dfprob_region(val,markers,region):
  # Calculate detecion probability for a region
  # alowed placements are where center of gene conversion tract is in [lbd,lbd[ of the region
  ExpVal=0

  # Find markers overlapping the region
  start=region.lbd-floor((val-1)/2) #included
  end  =region.ubd+floor((val)/2) #not included
  rmrks=markers[(markers.Chr==region.Chr) & (markers.pos>=start) & (markers.pos<end)]

  # Populate list according to first placement of gene conversion tract of size val
  # Find markers in list
  pre=rmrks[rmrks.pos<=start+val]
  thelist=list(pre.pos-start+1)

  # get diffs
  lmrks=rmrks[rmrks.pos>start+val]
  diffs=list(lmrks['diff'])

  # fix 1st diff
  diffs[0]=lmrks.iloc[0].pos-start-val

  # add one extra diff to stop at
  diffs.append(end-lmrks.iloc[-1].pos)
  
  # Slide a window of size val over region to calculate sum of expected #found gene conversion tracts
  # over all placements the gene conversion tract
  for d in list(rmrks[rmrks.pos>=region.lbd]['diff']):
    # The window is empty skip ahead until next marker enters window 
    # this case can only occur for first marker in region
    if not thelist:
      thelist.append(val)
      continue
    minv=min(thelist)

    # Check for markers exiting window before next marker enters it
    while minv<d:
      ExpVal+=minv*(1-(1-penetrance)**len(thelist))
      thelist.remove(minv)
      if not thelist:
        break
      thelist=[v-minv for v in thelist]
      d-=minv
      minv=min(thelist)
    # Skip ahead until next marker enters window
    ExpVal+=d*(1-(1-penetrance)**len(thelist))
    thelist=[v-d for v in thelist]
    if 0 in thelist: # next marker and enters window as the first marker in window exits
      thelist.remove(0)
    thelist.append(val)

  # We don't need to remove remaining elements as the last diff signals when to stop
#   while thelist:
#     minv=min(thelist)
#     ExpVal+=minv*(1-(1-penetrance)**len(thelist))
# 
#     thelist.remove(minv)
#     thelist=[v-minv for v in thelist]

  return ExpVal

def calc_detfun(diffs):
  plist=[dfprob(val,diffs) for val in klist]

  return plist

def calc_detfun_from_regions(markers,regions):
  plist=[sum([dfprob_region(val,markers,reg) for reg in regions.itertuples()]) for val in klist]

  return plist


parser = argparse.ArgumentParser(description='Calculate detection function from informative markers writing the output to [markerfile].p[PENETRANCE].detdf - The ouput file is suitable as input file for EMmodel.py')
parser.add_argument('markerfile', help='informative marker file')  
parser.add_argument('--regions',type=str,help="file with regions where the detecion functions should be calculated")
parser.add_argument('--dense',action='store_true',help="use denser grid")
parser.add_argument('--penetrance', type=float, default=0.5, help='probability of NCO marker being gene converted')  
parser.add_argument('--output', type=str, help='name of output file')
args=parser.parse_args()

penetrance=args.penetrance
thefile=args.markerfile

print('Processing:', thefile)

markers=read_markerfile(thefile)
genome=flatten_genome(markers)


if thefile.endswith('.gz'):
    thefile=thefile[:-3]

if args.output:
    detdfpath=args.output
else:
    detdfpath=thefile+".p%s.detdf"%penetrance

if args.dense:
    klist=list(range(1,10,1))
    klist += list(range(10,100,5))
    klist += list(range(100,500,10))
    klist += list(range(500,1000,50))
    klist += list(range(1000,10000,500))
    klist += list(range(10000,100001,5000))
    detdfpath=thefile+".p%s.dense.detdf"%penetrance

print('writing output to',detdfpath)

if args.regions:
    regions=pd.read_csv(args.regions,sep='\t')
    plist=calc_detfun_from_regions(markers,regions)
else:
    plist=calc_detfun(list(markers['diff']))

detdf = pd.DataFrame({"k":klist,"p":plist})
detdf.to_csv(detdfpath,sep="\t",index=False)




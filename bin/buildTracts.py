#!/usr/bin/python3
import pandas as pd
import sys
import argparse

def printFormatted(tr,mrk):
    for row in tr.itertuples():
        chrom=row.Chromosome
        lbd=row.First_converted_SNP
        ubd=row.Last_converted_SNP

        lbfilters=[mrk.Chr==chrom, mrk.pos>=lbd-100000, mrk.pos<lbd] 
        lgap=list(reversed(list(lbd-mrk[lbfilters[0] & lbfilters[1] & lbfilters[2]]['pos'])))
        lgstr=','.join(str(s) for s in lgap[:20])

        ubfilters=[mrk.Chr==chrom, mrk.pos<=ubd+100000, mrk.pos>ubd] 
        rgap=list(mrk[ubfilters[0] & ubfilters[1] & ubfilters[2]]['pos']-ubd)
        rgstr=','.join(str(s) for s in rgap[:20])

        print(chrom,lbd,ubd,lgstr,rgstr,sep="\t")

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='builds tracts file from GeneConversion csv')
    parser.add_argument('csv', type=str, help='file containing the tracts')
    parser.add_argument('--markerdir', type=str, help='directory containing marker files')
    args = parser.parse_args()

    tracts=pd.read_csv(args.csv,delimiter=',')

    t2=tracts[tracts['Generation']=="F2"]
    t5=tracts[tracts['Generation']!="F2"]

    print("Chr","lbd","ubd","lgap","rgap",sep='\t')

    F2=pd.read_csv('%s/F2.markerset'%args.markerdir,delimiter='\t')
    printFormatted(t2,F2)

    F45=pd.read_csv('%s/F45.markerset'%args.markerdir,delimiter='\t')
    printFormatted(t5,F45)


import pandas as pd
import argparse

label2markers={}
def loadMarkers(label):
    label2markers[label]=pd.read_csv('%s/%s.markerset'%(args.markerdir,label),delimiter='\t')


def main(df):
    print("Chromosome","First_converted_SNP","Last_converted_SNP","All_converted_SNPs","Number_of_converted_SNPs","Complex_NCO_Flag","Complex_CO_Flag","gcBetween","allBetween",sep='\t')
    for i in df.itertuples():
        label=i.Generation
        if i.Generation!='F2':
            label="F45"
        if label not in label2markers:
            loadMarkers(label)
        mrk=label2markers[label]
        between=sum((mrk['Chr']==i.Chromosome) & (mrk['pos']>i.First_converted_SNP) & (mrk['pos']<i.Last_converted_SNP))
        gcbetween=max(0,i.Number_of_converted_SNPs-2)
        print(i.Chromosome,i.First_converted_SNP,i.Last_converted_SNP,i.All_converted_SNPs,i.Number_of_converted_SNPs,i.Complex_NCO_Flag,i.Complex_CO_Flag,gcbetween,between,sep='\t')



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Find non converted markers within gene conversions')
    parser.add_argument('--markerdir', type=str, help='directory containing marker files')
    parser.add_argument('csv', type=str, help='file containing the tracts')
    args = parser.parse_args()

    df=pd.read_csv(args.csv,delimiter=',')
    main(df)


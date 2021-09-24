# coding: utf-8
import pandas as pd
import sys
import os
predSVFile=sys.argv[1]
outputFile=sys.argv[2]
fromChromosomeFormat=sys.argv[3] #either chrX or X, chr1, or 1, ..etc
print("input arguments were (three are required): ",sys.argv)
df = pd.read_csv(predSVFile, sep="\t")
df.head()
df.columns
rename={'#chrA':'ChrA','posA':'PosA', 'ortA':'OrtA','chrB':'ChrB', 'posB':'PosB', 'ortB':'OrtB'}
df = df.rename(columns=rename)

#### make sure no non-numeric values
ix1 = pd.to_numeric(df['PosA'], errors='coerce').isnull()
df = df.loc[~ix1]
assert(df.loc[ix1].shape[0]==0),"non-numeric values detected in the swab file in PosA"

ix2 = pd.to_numeric(df['PosB'], errors='coerce').isnull()
df = df.loc[~ix2]
assert(df.loc[ix2].shape[0]==0),"non-numeric values detected in the swab file in PosB"


## change the name of chromosomes from chr1 to 1 
def addChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==2 else 'chr'+row.strip()
    else:
        new_r = str(int(row)).strip()
        return 'chr'+new_r

def removeChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==1 else row.split("chr")[1].strip()
    else:
        new_r = str(int(row)).strip()
        return new_r
def formatChrom(row, fromChromosomeFormat):
    if len(fromChromosomeFormat.split("chr"))>1: #from format is chrX, ..etc
        return addChr(row)
    else:
        return removeChr(row)
        
df.loc[:,'ChrA'] = df['ChrA'].apply(lambda r: formatChrom(r,fromChromosomeFormat)).astype(pd.StringDtype())
df.loc[:,'ChrB'] = df['ChrB'].apply(lambda r: formatChrom(r,fromChromosomeFormat)).astype(pd.StringDtype())

#### if from format is X (i.e. we are lifting over from hg19, chromosome names such as Un_g1000224, ..etc are not acceptable and will cause liftover_sv.sh to fail
if len(fromChromosomeFormat.split("chr")) == 1: #lifting from hg19
    ## detecting when chromosome name is not X or Y and is not numeric
    idx1 = (~df['ChrA'].isin(['X','Y']) ) & (~df['ChrA'].str.isnumeric()) 
    df = df.loc[~idx1]
    idx2 = (~df['ChrB'].isin(['X','Y']) ) & (~df['ChrB'].str.isnumeric()) 
    df = df.loc[~idx2]

##### make sure that positions are intgers
df.loc[:,'PosA'] = df['PosA'].astype(int)
df.loc[:,'PosB'] = df['PosB'].astype(int)

df.head()
df.to_csv(outputFile, sep="\t", index=False)
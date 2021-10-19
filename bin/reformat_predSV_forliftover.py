# coding: utf-8
import pandas as pd
import sys
import os
predSVFile=sys.argv[1]
outputFile=sys.argv[2]
fromChromosomeFormat=sys.argv[3] #either chrX or X, chr1, or 1, ..etc
inputFormat=sys.argv[4] #either rawCrest or bedpe
print("input arguments were (four are required): ",sys.argv)
if inputFormat == "rawCrest":
    df = pd.read_csv(predSVFile, sep="\t")
    rename={'#chrA':'ChrA','posA':'PosA', 'ortA':'OrtA','chrB':'ChrB', 'posB':'PosB', 'ortB':'OrtB'}
    df = df.rename(columns=rename)
elif inputFormat == "bedpe":
	df = pd.read_csv(predSVFile, sep="\t", header=None)
	### from bedpe to rawCrest format one needs to get the coords 1-based
	### here take the padded coord instead
	df = df.iloc[:,[0,2,3,5,6,7,8,9]]
	cols = ['ChrA','PosA','ChrB','PosB','Type','score','OrtA','OrtB']
	df.columns = cols
	df = df[['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','score']]
	

#### drop null lines
idx = (~df['ChrA'].isnull()) & (~df['ChrB'].isnull()) & (~df['PosA'].isnull()) & (~df['PosB'].isnull())
df = df.loc[idx,:]

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

#### if from format is chrX (i.e. we are lifting over from hg38, chromosome names with alt in them (e.g. chr7_KI270803v1_alt will cause liftover_sv.sh to fail b/c these are not in the reference /research/rgs01/resgen/ref/tartan/runs/ad_hoc/NAME-4tmOgtgc/output/sequence/fasta/GRCh38_no_alt.fa
if len(fromChromosomeFormat.split("chr")) > 1: #lifting from hg38
    ## detecting when chromosoe name has _alt in it
    idx1 = df['ChrA'].apply(lambda r: True if len(r.split("_alt"))>1 else False)
    df = df.loc[~idx1]
    idx2 = df['ChrB'].apply(lambda r: True if len(r.split("_alt"))>1 else False)
    df = df.loc[~idx2]

##### make sure that positions are intgers
df.loc[:,'PosA'] = df['PosA'].astype(int)
df.loc[:,'PosB'] = df['PosB'].astype(int)

df.head()
df.to_csv(outputFile, sep="\t", index=False)

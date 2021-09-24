#!/usr/bin/python3
import numpy as np
import pandas as pd
import os
import sys
import time

"""
    takes the file  "$p"_hg38lo.swap.tab which is output of the liftover_flatfile.pl and prepares it by:
    1. converting positions to numeric values
    2. updating end positions to start+1
"""

inputFWpath = sys.argv[1]
outputFWpath = sys.argv[2]

df = pd.read_csv(inputFWpath, sep='\t')
print(df.head())
assert(df.columns.isin('ChrA\tPosA\tOrtA\tSClipA\tChrB\tPosB\tOrtB\tSClipB\tType\tCoverA\tCoverB\tConsLenA\tConsLenB\tMeanPIDA\tPctRepA\tMeanPIDB\tPctRepB\tConsMapStart\tConsMapStartChr\tConsMapStartPos\tConsMapEnd\tConsMapEndChr\tConsMapEndPos\tConsensusSeq\tInsSeq\tMicroHomo\tGeneA\tGeneB\tInvConsMapStart\tInvConsMapStartChr\tInvConsMapStartPos\tInvConsMapEnd\tInvConsMapEndChr\tInvConsMapEndPos\tInvConsensusSeq\tInvInsSeq\tInvMicroHomo\tliftover_ok_PosA\tliftover_chr_PosA\tliftover_base_PosA\tliftover_ort_PosA\tliftover_ok_PosB\tliftover_chr_PosB\tliftover_base_PosB\tliftover_ort_PosB\tliftover_min_match\tliftover_interval_length_delta\tliftover_interval_length_delta_fraction'.split("\t")).sum()==48),"some columns are missing! check the dataframe columns{}".format(df.columns.values)

##### take only the liftover columns and assign them simple names
dflo = df[['liftover_chr_PosA','liftover_base_PosA','liftover_ort_PosA','liftover_chr_PosB','liftover_base_PosB','liftover_ort_PosB','Type','CoverA','CoverB']]
dflo.columns = ['chrom1','start1','strand1','chrom2','start2','strand2','name','cover1','cover2']

#### make sure no non-numeric values
ix1 = pd.to_numeric(dflo['start1'], errors='coerce').isnull()
dflo = dflo.loc[~ix1]
assert(dflo.loc[ix1].shape[0]==0),"non-numeric values detected in the swab file in start1"

ix2 = pd.to_numeric(dflo['start2'], errors='coerce').isnull()
dflo = dflo.loc[~ix2]
assert(dflo.loc[ix2].shape[0]==0),"non-numeric values detected in the swab file in start2"

##### convert positions to integers
dflo.loc[:,'start1'] = dflo['start1'].apply(lambda r: str(int(r)).strip()).astype(pd.StringDtype())
dflo.loc[:,'start2'] = dflo['start2'].apply(lambda r: str(int(r)).strip()).astype(pd.StringDtype())
dflo = dflo.astype({'start1':pd.Int64Dtype(),'start2':pd.Int64Dtype()})

### change the name of chromosomes from 1 to chr1 
def formatChrom(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==2 else 'chr'+row.strip()
    else:
        new_r = str(int(row)).strip()
        return 'chr'+new_r
        
dflo.loc[:,'chrom1'] = dflo['chrom1'].apply(lambda r: formatChrom(r)).astype(pd.StringDtype())
dflo.loc[:,'chrom2'] = dflo['chrom2'].apply(lambda r: formatChrom(r)).astype(pd.StringDtype())
#dflo.loc[:,'chrom1'] = dflo['chrom1'].apply(lambda r: str(r).strip() if len(str(r).split('chr'))==2 else 'chr'+str(int(r)).strip()).astype(pd.StringDtype())
#dflo.loc[:,'chrom2'] = dflo['chrom2'].apply(lambda r: str(r).strip() if len(str(r).split('chr'))==2 else 'chr'+str(int(r)).strip()).astype(pd.StringDtype())

#### update ends so that it is ready for pairToPair
dflo.loc[:,'end1'] = dflo['start1'].values + 1
dflo.loc[:,'end2'] = dflo['start2'].values + 1


##### add a summary coverage and reformat to a bedpe file
dflo.loc[:,'geomMeanCov'] = np.sqrt(dflo['cover1']*dflo['cover2']).values

dflo = dflo['chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tcover1\tcover2'.split('\t')]
assert(np.all(dflo.columns == 'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tcover1\tcover2'.split('\t'))), "subsetting columns is not successfull"

#### copy over the first step
dflo.to_csv(outputFWpath, sep='\t', header=None, index=False)

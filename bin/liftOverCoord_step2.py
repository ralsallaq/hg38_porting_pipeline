#!/usr/bin/python3
import numpy as np
import pandas as pd
import os
import sys
import time

inputFWpath = sys.argv[1]
outputFWpath = sys.argv[2]

df = pd.read_csv(inputFWpath, sep='\t')
print(df.head())
assert(df.columns.isin('chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tInvstart1\tInvstart2'.split('\t'))
                                         .sum()==12),"some columns are missing! check the dataframe columns {}".format(df.columns.values)

##### convert positions to integers
df = df.astype({'start1':pd.Int64Dtype(),'end1':pd.Int64Dtype(),'start2':pd.Int64Dtype(),
                'end2':pd.Int64Dtype(),'Invstart1':pd.Int64Dtype(),'Invstart2':pd.Int64Dtype()})

#### make sure no non-numeric values
ix1 = pd.to_numeric(df['start1'], errors='coerce').isnull()
assert(df.loc[ix1].shape[0]==0),"non-numeric values detected in the swab file in start1"

ix2 = pd.to_numeric(df['start2'], errors='coerce').isnull()
assert(df.loc[ix2].shape[0]==0),"non-numeric values detected in the swab file in start2"

#### update ends so that it is ready for pairToPair
df.loc[:,'Invend1'] = df['Invstart1'].values + 1
df.loc[:,'Invend2'] = df['Invstart2'].values + 1
df = df['chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tInvstart1\tInvend1\tInvstart2\tInvend2'.split('\t')]

df = df['chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tInvstart1\tInvstart2'.split('\t')]
assert(np.all(df.columns == 'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomMeanCov\tstrand1\tstrand2\tInvstart1\tInvstart2'.split('\t'))), "subsetting columns is not successfull"

#### copy over the first step
df.to_csv(outputFWpath, sep='\t', index=False)

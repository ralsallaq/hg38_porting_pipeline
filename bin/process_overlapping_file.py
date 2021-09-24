#!/usr/bin/python3
import numpy as np
import pandas as pd
import os
import sys
import time


inputFWpath = sys.argv[1]
df = pd.read_csv(inputFWpath, sep='\t', header=None)
#'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tcovAcovB\tstrand1\tstrand2\tInvstart1\tInvstart2'.split('\t')
#print(df.head())
df.loc[:,'match-type'] = (df.loc[:,6]==df.loc[:,16]).astype(int)
df.loc[:,'absDiffcovAcovB'] = abs(df.loc[:,7]-df.loc[:,17]).values
#print(df.head())
print("overlap-match-type","\toverlap-discordant-type","\tmeanAbsDiffcovAcovB")
print(df.loc[:,'match-type'].sum(),"\t",df.shape[0]-df.loc[:,'match-type'].sum(),"\t",df.loc[:,'absDiffcovAcovB'].mean())


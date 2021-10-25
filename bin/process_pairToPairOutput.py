#!/usr/bin/python3
import numpy as np
import pandas as pd
import os
import sys
import time


inputFWpath = sys.argv[1]
df = pd.read_csv(inputFWpath, sep='\t', header=None)
#'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tgeomCov\tstrand1\tstrand2\tInvstart1\tInvstart2'.split('\t')
#print(df.head())
if df.columns.shape[0]==20: #dbedpe-->overlapping results
    df.loc[:,'match-type'] = (df.loc[:,6]==df.loc[:,16]).astype(int)
    df.loc[:,'absDiffGeomCov'] = abs(df.loc[:,7]-df.loc[:,17]).values
    df.loc[:,'midGeomCov'] = (df.loc[:,7]+df.loc[:,17]).values/2
    #print(df.head())
    print("nu_overalps","\toverlap-match-type","\toverlap-discordant-type","\tmeanAbsDiffgeomCov","\tmeanMidGeomCov")
    print(df.shape[0],"\t", df.loc[:,'match-type'].sum(),"\t",df.shape[0]-df.loc[:,'match-type'].sum(),"\t",df.loc[:,'absDiffGeomCov'].mean(),"\t",df.loc[:,'midGeomCov'].mean())
elif df.columns.shape[0]<20 and df.columns.shape[0]>=10: #bedpe--> missing/extra results
    df.loc[:,'geomCov'] = df.loc[:,7]
    #print(df.head())
    print("meangeomCov")
    print(df.loc[:,'geomCov'].mean())
else:
    print("not valid input file, input should be either dbedpe (overlapping results) or bedpe (missing/extra results)")
    sys.exit(1)

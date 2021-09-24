/*============================= miscellaneous utility functions  ===============================*/
// to convert predSV raw files to bedpe
process reformatLOPredSVToBEDPE {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_${params.FromTo}", mode: 'copy'

    input:
    tuple val(anlsRunName), val(pairName), val(genome), path(rawSomaticSVcalls)

    output:
    tuple val(anlsRunName), val(pairName), val(genome), path("${pairName}_${genome}.bedpe"), emit: bedpeCrestOutput_ch

    beforeScript "export PATH=${workflow.projectDir}/bin:$PATH; source set_env.sh"

    """
#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys, os
import glob
#        takes a text file *.predSV which is output of liftover_sv.sh
#        and reformat it to bedpe file 
print('reading in predSV file ','${rawSomaticSVcalls}')
df = pd.read_csv('${rawSomaticSVcalls}', sep='\\t')
#### geometric mean of the coverage of the two points
df.loc[:,'geomMeanCov'] = np.sqrt(df['CoverA']*df['CoverB']).values
#ChrA    PosA    OrtA    ChrB    PosB    OrtB    Type    CoverA  CoverB
cols={'ChrA':'chrom1','ChrB':'chrom2','OrtA':'strand1','OrtB':'strand2','PosA':'start1','PosB':'start2','Type':'name','CoverA':'cover1','CoverB':'cover2'}
df.rename(columns=cols, inplace=True)

#### convert position to integer
### first drop null values 
df =df[~df['start1'].isnull()]
df =df[~df['start2'].isnull()]

#### what is left are either numbers (float/int) or strings
#### convert to str
df.loc[:,'start1'] = df['start1'].astype(str)
df.loc[:,'start2'] = df['start2'].astype(str)

def isfloat(x):
    try:
        a = float(x)
    except (TypeError, ValueError):
        return False
    else:
        return True

def isint(x):
    try:
        a = float(x)
        b = int(a)
    except (TypeError, ValueError):
        return False
    else:
        return a == b

### drop non_numeric (not int not float)
idx1 = df['start1'].apply(lambda r: True if isfloat(r) or isint(r) else False) 
df = df.loc[idx1]
idx2 = df['start2'].apply(lambda r: True if isfloat(r) or isint(r) else False) 
df = df.loc[idx2]

df.loc[:,'start1'] = df['start1'].astype(float).astype(int)
df.loc[:,'start2'] = df['start2'].astype(float).astype(int)

df.loc[:,'end1'] = df['start1']+1
df.loc[:,'end2'] = df['start2']+1
df = df[['chrom1', 'start1','end1','chrom2', 'start2','end2','name', 'geomMeanCov','strand1','strand2','cover1', 'cover2']]


### change the name of chromosomes from 1 to chr1 
def formatChrom(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==2 else 'chr'+row.strip()
    else:
        new_r = str(int(float(row))).strip()
        return 'chr'+new_r
        
df.loc[:,'chrom1'] = df['chrom1'].apply(lambda r: formatChrom(r)).astype(pd.StringDtype())
df.loc[:,'chrom2'] = df['chrom2'].apply(lambda r: formatChrom(r)).astype(pd.StringDtype())

df.to_csv('${pairName}_${genome}.bedpe',sep='\\t', header=None, index=False)

print('done')

    """
  }

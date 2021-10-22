/*============================= miscellaneous utility functions  ===============================*/
// to convert predSV raw files to bedpe
process reformatLOPredSVToBEDPE {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_${params.FromTo}", mode: 'copy'

    input:
    tuple val(anlsRunName), val(pairName), val(genome), path(rawSomaticSVcalls)
    val(finalGenome)
    
    output:
    tuple val(anlsRunName), val(pairName), val(genome), path("${pairName}_${genome}_In${finalGenome}.bedpe"), emit: bedpeCrestOutput_ch

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
if df.columns.isin(['CoverA','CoverB']).sum()==2:
    df.loc[:,'score'] = np.sqrt(df['CoverA']*df['CoverB']).values
elif df.columns.isin(['score']).sum()==1:
    pass
    
cols={'ChrA':'chrom1','ChrB':'chrom2','OrtA':'strand1','OrtB':'strand2','PosA':'start1','PosB':'start2','Type':'name'}
df.rename(columns=cols, inplace=True)

#### convert position to integer
### first drop null values 
idx = (~df['chrom1'].isnull()) & (~df['chrom2'].isnull()) & (~df['start1'].isnull()) & (~df['start2'].isnull())
df = df.loc[idx,:]

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

#### because predSV are 1-based
df.loc[:,'start1'] = df['start1'].astype(float).astype(int) - 1
df.loc[:,'start2'] = df['start2'].astype(float).astype(int) - 1

df.loc[:,'end1'] = df['start1']+1
df.loc[:,'end2'] = df['start2']+1
df = df[['chrom1', 'start1','end1','chrom2', 'start2','end2','name', 'score','strand1','strand2']]


### change the name of chromosomes from 1 to chr1 or chr1 to 1 depending on the finalFormat
finalFormat = "chrX" if '${finalGenome}'=="hg38" else "X"

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

def formatChrom(row, finalFormat):
    if len(finalFormat.split("chr"))>1: #final format is chrX, ..etc
        return addChr(row)
    else:
        return removeChr(row)
        
df.loc[:,'chrom1'] = df['chrom1'].apply(lambda r: formatChrom(r,finalFormat)).astype(pd.StringDtype())
df.loc[:,'chrom2'] = df['chrom2'].apply(lambda r: formatChrom(r,finalFormat)).astype(pd.StringDtype())

df.to_csv('${pairName}_${genome}_In${finalGenome}.bedpe',sep='\\t', header=None, index=False)

print('done')

    """
  }

//to combine predSV files from the same genome for all sample pairs into one
process concatFusions {
    label 'io_mem'
    publishDir "${params.outD}/crest-post/", mode: 'copy'

    input:
    path(fusionSVFiles)
    val(genome)

    output:
    path("fusionSVFromAllpairs_${genome}.csv"), emit: combinedCrestPostOutput_ch

    """
#!/usr/bin/env python3
import pandas as pd
fusionFiles = "${fusionSVFiles}".split(" ")
pairs = [f.split("-")[0].strip() for f in fusionFiles]
df_files = pd.DataFrame({'file':fusionFiles,'pair':pairs})
#sort by pairs so that files are aligned regardless to genome
df_files = df_files.sort_values(by='pair')
print(df_files)

#read files into DFs
dfs = []
for i, row in df_files.iterrows(): 
    temp_df = pd.read_csv(row['file'], sep="\t")
    dfs.append(temp_df)

dfs_all = pd.concat(dfs, axis=0)

dfs_all.to_csv("fusionSVFromAllpairs_"+"${genome}"+".csv", index=False)
    """
}

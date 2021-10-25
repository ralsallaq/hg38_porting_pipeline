/*============================= miscellaneous utility functions  ===============================*/
// to convert predSV raw files to bedpe

process reformatForLiftOver{
    label 'io_mem'
    publishDir "${params.outD}/${params.OutputDir}/", mode: 'copy'

    input:
    //path(inputFile)
    tuple val(anlsRunName), val(pairName), val(genome), path(inputFile)
    val(inputFormat)
    val(fromGenome)

    output:
    //path("readyForLO_from_${fromGenome}.predSV"), emit: reformattedForLiftOver_ch
    tuple val(anlsRunName), val(pairName), val(genome), path("${pairName}_readyForLO_from_${fromGenome}.${inputFormat}"), emit: reformattedForLiftOver_ch 

    """
#!/usr/bin/env python3
import pandas as pd
import sys
import os
if '${fromGenome}'=='hg38':
    fromChromosomeFormat='chrX'
elif '${fromGenome}'=='hg19':
    fromChromosomeFormat='X'
else:
    print("fromGenome can be only hg38 or hg19 ..exiting!")
    sys.exit(1)

inputFormat = '${inputFormat}'
inputFile= '${inputFile}'

if inputFormat == "rawCrest":
    df = pd.read_csv(inputFile, sep="\\t")
    rename={'#chrA':'ChrA','posA':'PosA', 'ortA':'OrtA','chrB':'ChrB', 'posB':'PosB', 'ortB':'OrtB'}
    df = df.rename(columns=rename)
elif inputFormat == "bedpe":
    df = pd.read_csv(inputFile, sep="\\t", header=None)
    ### from bedpe to rawCrest format one needs to get the coords 1-based
    ### here take the padded coord instead
    df = df.iloc[:,[0,2,3,5,6,7,8,9]]
    cols = ['ChrA','PosA','ChrB','PosB','Type','score','OrtA','OrtB']
    df.columns = cols
    df = df[['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','score']]
elif inputFormat == "crestPost":
    df = pd.read_csv(inputFile, sep="\\t")
    rename={'OrientA':'OrtA','OrientB':'OrtB'}
    df = df.rename(columns=rename)
    if df.columns.isin(['rating', 'Usage']).sum()==2:
        cols=['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB','Fusion Gene','rating','Usage','Sample(s)']
    elif df.columns.isin(['rating']).sum()==1:
        cols=['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB','Fusion Gene','rating','Sample(s)']
    elif df.columns.isin(['Usage']).sum()==1:
        cols=['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB','Fusion Gene','Usage','Sample(s)']
    else:
        cols=['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB','Fusion Gene','Sample(s)']

    df = df[cols]
    if '${params.keepAnnotatedOnly}'=='true':
        #drop SVs without Fusion Genes
        df = df[~df['Fusion Gene'].isnull()]

	

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
#df.to_csv('readyForLO_from_${fromGenome}.predSV', sep="\\t", index=False)
df.to_csv('${pairName}_readyForLO_from_${fromGenome}.${inputFormat}', sep="\\t", index=False)

print('done')
    """
}

process reformatLOPredSVToBEDPE {
    label 'io_mem'
    publishDir "${params.outD}/${params.OutputDir}/", mode: 'copy'

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

bedpe_cols = ['chrom1', 'start1','end1','chrom2', 'start2','end2','name', 'score','strand1','strand2']

cols = bedpe_cols + df.columns[~df.columns.isin(bedpe_cols)].tolist()

#df = df[['chrom1', 'start1','end1','chrom2', 'start2','end2','name', 'score','strand1','strand2']]
df = df[cols]


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
    publishDir "${params.outD}/${params.OutputDir}/", mode: 'copy'

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
    temp_df = pd.read_csv(row['file'], sep="\\t")
    dfs.append(temp_df)

dfs_all = pd.concat(dfs, axis=0)

dfs_all.to_csv("fusionSVFromAllpairs_"+"${genome}"+".csv", index=False, sep="\\t")
    """
}

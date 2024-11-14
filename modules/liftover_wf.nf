/*============================= liftover workflow ===============================*/

//to liftover SVs
process liftoverSVs {
    label 'io_mem'
    publishDir "${params.outD}/predSV_${params.FromTo}", mode: 'copy'

    input:
    tuple val(anlsRunName), val(pairName), val(genome), path(rawSomaticSVcalls)
    val(loFromToOfficial)
    val(lofromto)
    val(lofromFromat)
    val(inputFormat)

    output:
    tuple val(anlsRunName), val(pairName), val(genome), path("${pairName}_${genome}_${lofromto}.${inputFormat}"), emit: liftoverBEDPE_ch
    
   // beforeScript "export PATH=${workflow.projectDir}/bin:$PATH; source set_env.sh"

    """
#!/usr/bin/bash

echo "you need to set up the environment for liftover_flatfile.pl and other programs by: \
cbload seq_anls_ops \
cbload snv-annovar \
module load vcftools/0.1.15 \
module load samtools/1.15.1 \
module load sjcb/ucsc/012314 \
module load perl/5.10.1 #sjcb perl 5.30.0 did not work \
1) vcftools/0.1.15   2) samtools/1.15.1   3) sjcb/ucsc/012314   4) perl/5.10.1"


echo "official names of genome from and to are needed for liftover_sv.sh"

lofromOfficial=`echo "${loFromToOfficial}" | cut -f1 -d "|"`
lotoOfficial=`echo "${loFromToOfficial}"|cut -f3 -d "|"`
lofrom=`echo "${lofromto}"|cut -f1 -d "_"`
loto=`echo "${lofromto}"|cut -f3 -d "_"`

echo "processing liftover from \$lofrom to \$loto"
echo "processing liftover from offcial \$lofromOfficial to \$lotoOfficial"

### prep file
echo "prep the raw SV calls file for liftover_sv.sh script"

reformattedFile="${lofromto}"_reformatted.predSV

python  ${workflow.projectDir}/bin/reformat_predSV_forliftover.py ${rawSomaticSVcalls} \$reformattedFile ${lofromFromat} ${inputFormat}

##### apply leftover 
echo "applying liftover_sv.sh which produces the same *.predSV file with extra columns"

liftover_sv.sh -w tempDir \$lofromOfficial \$lotoOfficial  \$reformattedFile ${pairName}_${genome}_${lofromto}_raw

wait
echo "save a predSV formatted file"
python -c " \

import numpy as np
import pandas as pd
import os
import sys
import time

inputFWpath = '${pairName}_${genome}_${lofromto}_raw'
outputFWpath = os.path.basename(inputFWpath).split('_raw')[0]+'.${inputFormat}'

dflo = pd.read_csv(inputFWpath, sep='\t')
print(dflo.head())
assert(dflo.columns.isin('ChrA\tPosA\tOrtA\tChrB\tPosB\tOrtB\tType\tliftover_ok_PosA\tliftover_chr_PosA\tliftover_base_PosA\tliftover_ort_PosA\tliftover_ok_PosB\tliftover_chr_PosB\tliftover_base_PosB\tliftover_ort_PosB\tliftover_min_match\tliftover_interval_length_delta\tliftover_interval_length_delta_fraction'.split('\t')).sum()==18),'some columns are missing! check the dataframe columns{}'.format(dflo.columns.values)

##### take only the liftover columns and assign them simple names
if '${inputFormat}'=='rawCrest':
    dflo = dflo[['liftover_chr_PosA','liftover_base_PosA','liftover_ort_PosA','liftover_chr_PosB','liftover_base_PosB','liftover_ort_PosB','Type','CoverA','CoverB']]
    dflo.columns = ['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB']
elif '${inputFormat}'=='bedpe':
    dflo = dflo[['liftover_chr_PosA','liftover_base_PosA','liftover_ort_PosA','liftover_chr_PosB','liftover_base_PosB','liftover_ort_PosB','Type','score']]
    dflo.columns = ['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','score']
elif '${inputFormat}'=='crestPost':
    if dflo.columns.isin(['rating', 'Usage']).sum()==2:
        addcols=['Fusion Gene','rating','Usage']
    elif dflo.columns.isin(['rating']).sum()==1:
        addcols=['Fusion Gene','rating']
    elif dflo.columns.isin(['Usage']).sum()==1:
        addcols=['Fusion Gene','Usage']
    else:
        addcols=['Fusion Gene']

    dflo = dflo[['liftover_chr_PosA','liftover_base_PosA','liftover_ort_PosA','liftover_chr_PosB','liftover_base_PosB','liftover_ort_PosB','Type','CoverA','CoverB'] + addcols]
    dflo.columns = ['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','CoverA','CoverB'] + addcols



print('saving a liftover file')
dflo.to_csv(outputFWpath, sep='\t', index=False)

"
### end of python script

echo "done"

    """
}

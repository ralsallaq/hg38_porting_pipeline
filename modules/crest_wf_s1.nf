/*============================= crest workflow ===============================*/
workflow crest_wf {
    take:crest_runs_ch

    main:

    //run reformatToBEDPE
    //crest_runs_ch.subscribe { println "value: $it"} // this seems not to consume the channel!!

    reformatToBEDPE(crest_runs_ch)

    if (params.lift_over_fromTo == 'hg19_to_hg38') {
        //pass only pedpe for hg19
        liftoverCoord(reformatToBEDPE.out.filter{it[2]=='hg19'})

        //combine channels hg38 first and hg19 second
        // e.g. SJE2A001_D_G, MI1j7VRI, hg38, /research/groups/zhanggrp/home/ralsalla/test_crest_newConfig_hg38/crest_compare_newConfigJuly20_2021_wComb/work/7e/ce0c367e4fc15d064099647cc4b42e/SJE2A001_D_G_hg38.bedpe, jtkZLnf1, hg19, /research/groups/zhanggrp/home/ralsalla/test_crest_newConfig_hg38/crest_compare_newConfigJuly20_2021_wComb/work/42/7f906fcd4b8df3651682b2aa5d549a/SJE2A001_D_G_hg19lo
        bedpe_combined_ch = reformatToBEDPE.out.filter{it[2]=='hg38'}.map{it->[it[1], it[0], it[2], it[3]]}.join(liftoverCoord.out.map{it->[it[1], it[0], it[2], it[3]]})
        //reformatToBEDPE.out.filter{it[2]=='hg38'}.map{it->[it[1], it[0], it[2], it[3]]}.join(liftoverCoord.out.map{it->[it[1], it[0], it[2], it[3]]}).subscribe{println it}

        
        emit:
            reformatToBEDPE.out.filter{it[2]=='hg38'}.join(liftoverCoord.out)
    } else if (params.lift_over_fromTo == 'hg38_to_hg19') {
        //pass only pedpe for hg38
        liftoverCoord(reformatToBEDPE.out.filter{it[2]=='hg38'})

        //combine channels hg38 first and hg19 second
        bedpe_combined_ch = liftoverCoord.out.map{it->[it[1], it[0], it[2], it[3]]}.join(reformatToBEDPE.out.filter{it[2]=='hg19'}.map{it->[it[1], it[0], it[2], it[3]]})

        emit:
            liftoverCoord.out.join(reformatToBEDPE.out.filter{it[2]=='hg19'}).view()
    } else {
        log.error("Not valid lift_over_fromTo. Valid modes are: 'hg19_to_hg38' or 'hg38_to_hg19' only ")
        exit 1
    }


    //pass pairs of bedpe (one liftover and one not) to compare the two using pairToPair  
    compareSVs(bedpe_combined_ch)

    //summarize 
    //summarizeOverlapResults(compareSVs.out)
    compareSVs.out.map{it->it[1]}
        .collectFile(name: "${params.outD}/bedpe_comp_${params.lift_over_fromTo}/summary_all.out", keepHeader:true, skip:1, newLine:true)
        .view {file -> "pair\tNbreakpnts_hgFrom\tNbreakpnts_hgTo\tNoverlappingSVs\tNoverlap_match_type\tNoverlap_discordant_type\tNmissing\tNextra:${file.text}" }

    //get variant discordant rates for each pair for each sample used formula from functional equivalence paper (Nature Communincation 2018) 
    //formula for pair = (discordant + 0-only + 1-only + discordant_discordant_type)/(match+ discordant + match_discordant_type + discordant_discordant_type + 0-only + 1-only) 
    //match=overlapping positions AND same_strand AND sv-type is the same AND genotype is the same. match_discordant_type= all the same except sv-type. discordant=If there is overlap but genotypes are different irrespective if the sv-types are the same. discordant_discordant_type=if there is overlap and both genotype and sv-type are different.   
    
    
}//end of workflow

process reformatToBEDPE {
    label 'io_mem'
    publishDir "${params.outD}/bedpe", mode: 'copy'

    input:
    tuple val(anlsRunName), val(pairName), val(genome), file(rawSomaticSVcalls)

    output:
    tuple val(anlsRunName), val(pairName), val(genome), file("${pairName}_${genome}.bedpe"), emit: bedpeCrestOutput_ch

    """
#!/usr/bin/env python3
import pandas as pd
import sys, os
import glob
#        takes a text file *.predSV which is output of CREST
#        and reformat it to bedpe file 
print('reading in predSV file ','${rawSomaticSVcalls}')
df = pd.read_csv('${rawSomaticSVcalls}', sep='\\t')
df = df[['#chrA','ConsMapStartPos','chrB','ConsMapEndPos','Type','MeanPIDA','ortA','ortB','InvConsMapStartPos','InvConsMapEndPos']]
cols = cols={'#chrA':'chrom1','chrB':'chrom2','ortA':'strand1','ortB':'strand2','ConsMapStartPos':'start1','ConsMapEndPos':'start2','Type':'name','InvConsMapStartPos':'Invstart1','InvConsMapEndPos':'Invstart2'}
df.rename(columns=cols, inplace=True)
df.loc[:,'end1'] = df['start1']+1
df.loc[:,'end2'] = df['start2']+1
df = df[['chrom1', 'start1','end1','chrom2', 'start2','end2','name', 'MeanPIDA','strand1','strand2','Invstart1', 'Invstart2']]
df.loc[:,'chrom1'] = df.loc[:,'chrom1'].astype(str)
df.loc[:,'chrom2'] = df.loc[:,'chrom2'].astype(str)
df.loc[:,'chrom1'] = df['chrom1'].apply(lambda r: r if len(r.split('chr'))==2 else 'chr'+str(r))
df.loc[:,'chrom2'] = df['chrom2'].apply(lambda r: r if len(r.split('chr'))==2 else 'chr'+str(r))
df.to_csv('${pairName}_${genome}.bedpe',sep='\\t', header=None, index=False)

print('done')

    """
  }


process liftoverCoord {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_${params.lift_over_fromTo}", mode: 'copy'

    input:
    tuple val(anlsRunName), val(pairName), val(genome), file(rawSomaticSVcallsBEDPE)

    output:
    tuple val(anlsRunName), val(pairName), val(genome), file("${pairName}_${genome}lo_*.bedpe"), emit: bedpeliftover_ch
    
    //beforeScript "export PATH=${workflow.projectDir}/bin:$PATH; source set_env.sh"

    """
#!/usr/bin/bash
echo "you need to set up the environment for liftover_flatfile.pl and other programs by: \
cbload configs \
cbload common-scripts-internal \
"
module load ucsc/031411-legacy
module load python/3.7.0

lofrom=`echo "${params.lift_over_fromTo}"|cut -f1 -d "_"`
loto=`echo "${params.lift_over_fromTo}"|cut -f3 -d "_"`
echo "processing liftover from \$lofrom to \$loto"

if [ "${genome}" == "\$lofrom" ]; then
    ##### applying leftover 
    cat ${rawSomaticSVcallsBEDPE} | awk 'BEGIN{print "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tMeanPIDA\tstrand1\tstrand2\tInvstart1\tInvstart2"}1' >temp1 
    liftover_flatfile.pl -file temp1 -out ${pairName}_${genome}lo -field-ch chrom1 -field-ch chrom2 -field-pos start1 -field-pos start2 -genome-from ${genome} -genome-to \$loto -interval-interbase 

    wait
    echo "update the end positions to start+1 in *.swap.tab files and make sure the positions are all numeric"
    python ${workflow.projectDir}/bin/liftOverCoord_step1.py ${pairName}_${genome}lo.swap.tab ${pairName}_${genome}lo_.bedpe ${pairName}_invtemp 

    wait


    if [ -f ${pairName}_invtemp ]; then
        echo "INV file is detected will run liftover of the INV positions"
        liftover_flatfile.pl -file ${pairName}_invtemp -out ${pairName}_invs -field-ch chrom1 -field-ch chrom2 -field-pos Invstart1 -field-pos Invstart2 -genome-from ${genome} -genome-to \$loto -interval-interbase 
        wait
        #### now replace inv positions with liftover ones and update ends 
        python ${workflow.projectDir}/bin/liftOverCoord_step2.py ${pairName}_invs.swap.tab ${pairName}_${genome}lo_invs.bedpe
    
    fi

else 
    echo "assertion error! ${genome} != \$lofrom"
fi

    """
}


process compareSVs {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_comp_${params.lift_over_fromTo}", mode: 'copy'

    input:
    tuple val(pairName), val(anlsRunNameTo), val(genomeTo), file(bedpeFileTo), val(anlsRunNameFrom), val(genomeFrom), file(bedpeFileFrom)

    output:
    tuple val(pairName), file("summary_${pairName}_${params.lift_over_fromTo}"), file("${pairName}_wliftover_both_samestrand.dbedpe"), file("${pairName}_wliftover_missing.bedpe"), file("${pairName}_wliftover_extra.bedpe"), emit: pairToPairOutput_ch

    """
#!/usr/bin/bash

module load bedtools/2.30.0
module load python/3.7.0

printf "pair\tNbreakpnts_hgFrom\tNbreakpnts_hgTo\tNoverlappingSVs\tNoverlapping_match_type\tNoverlapping_discordant_type\tNmissing\tNextra\n" > summary_${pairName}_${params.lift_over_fromTo}
echo "the From genome is the current and the To genome is the new genome"

##### only the first 10 columns matter
echo "the bedpeFileFrom can be two files depending on whether there was _inv file or not in the liftover step"
if [ `echo "${bedpeFileFrom}"|cut -f 1 -d " "|grep -v '_invs'|wc -l` == 1 ]; then
    targetBEDPE=`echo "${bedpeFileFrom}"|cut -f 1 -d " "`
else
    targetBEDPE=`echo "${bedpeFileFrom}"|cut -f 2 -d " "`
fi

echo "target BEDPE file is set to \$targetBEDPE"

cat ${bedpeFileTo} |cut -f 1-10 > hgToF
echo "taking the From genome bedpe file with positions of breakpoint liftover to the To genome"
cat \$targetBEDPE|sed 1d |cut -f 1-10 > hgFromF
npnts_hgTo=`cat hgToF|wc -l`
npnts_hgFrom=`cat hgFromF|wc -l`

#### enforcing strandedness and existence in both

echo "a is From and b is To finding overlapping features in both"
echo "will save overlapping features for pair to ${pairName}_wliftover_both_samestrand.dbedpe which is a double bedpe (20 columns)"
pairToPair -f 1e-9 -type both  -is -a hgFromF -b hgToF > ${pairName}_wliftover_both_samestrand.dbedpe

python ${workflow.projectDir}/bin/process_overlapping_file.py ${pairName}_wliftover_both_samestrand.dbedpe >Nmatch_discordant_types

Noverlap_match_type=`cat Nmatch_discordant_types| sed 1d|cut -f 1`
Noverlap_discordant_type=`cat Nmatch_discordant_types| sed 1d|cut -f 2`

wait

## convert overlapping result to bedpe
cat ${pairName}_wliftover_both_samestrand.dbedpe | cut -f 1-10 > ${pairName}_wliftover_both_samestrand.bedpe

wait

echo "finding SVs found using From genome but ends do not overlap the ends of previous overlapping result with the To genome (i.e. in From but exclude the overlap with the To)"
pairToPair -f 1e-9 -type neither -a hgFromF -b ${pairName}_wliftover_both_samestrand.bedpe > ${pairName}_wliftover_missing.bedpe

wait

echo "finding SVs found using To genome but ends do not overlap the ends of previous overlapping result (i.e. in To but exclude the overlap with the From)"
pairToPair -f 1e-9 -type neither -a hgToF -b ${pairName}_wliftover_both_samestrand.bedpe  > ${pairName}_wliftover_extra.bedpe


wait
nu_missing=`cat ${pairName}_wliftover_missing.bedpe|wc -l`
nu_extra=`cat ${pairName}_wliftover_extra.bedpe|wc -l`

nu_overlaps=`cat ${pairName}_wliftover_both_samestrand.bedpe|wc -l`
printf "${pairName}\\t\$npnts_hgFrom\\t\$npnts_hgTo\\t\$nu_overlaps\\t\$Noverlap_match_type\\t\$Noverlap_discordant_type\\t\$nu_missing\\t\$nu_extra" >>  summary_${pairName}_${params.lift_over_fromTo}

echo "done"

    """
    
}

/*
process summarizeOverlapResults {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_comp_${params.lift_over_fromTo}", mode: 'copy'

    input:
    tuple val(pairName), file(summaryForPair), file(overlapForPair), file(missingForPair), file(extraForPair)

    output:
    file(summary_${params.lift_over_fromTo})

    """
#!/usr/bin/bash

    """
}
*/

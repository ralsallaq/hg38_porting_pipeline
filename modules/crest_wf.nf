/*============================= crest workflow ===============================*/
include { liftoverSVs as liftoverSVs_round1 }  from './liftover_wf'
include { liftoverSVs as liftoverSVs_round2 }  from './liftover_wf'
include { liftoverSVs }  from './liftover_wf'
include {reformatLOPredSVToBEDPE as refor_hg38_ToBEDPE } from './miscellaneous.nf'
include {reformatLOPredSVToBEDPE as refor_hg19_ToBEDPE } from './miscellaneous.nf'

workflow crest_wf {
    take:crest_runs_ch

    main:


    if ( params.FromTo == 'hg19_to_hg38') {

        /***** to reduce liftover artifacts liftover for hg38 data will be done twice  and for hg19 only once****/
        //liftover hg38 to hg19
        liftoverSVs_round1(crest_runs_ch.filter{it[2]=='hg38'}, channel.value('GRCh38|to|GRCh37-lite'), channel.value('hg38_to_hg19'), channel.value('chrX'))
        //liftover hg19 back to hg38
        liftoverSVs_round2(liftoverSVs_round1.out, channel.value('GRCh37-lite|to|GRCh38'), channel.value('hg19_to_hg38'), channel.value('X'))
        
        //reformat to bedpe
        refor_hg38_ToBEDPE(liftoverSVs_round2.out)

        //pass only predSV for hg19
        liftoverSVs(crest_runs_ch.filter{it[2]=='hg19'}, channel.value('GRCh37-lite|to|GRCh38'), channel.value('hg19_to_hg38'),channel.value('X'))

        //reformat to bedpe
        refor_hg19_ToBEDPE(liftoverSVs.out)

        //combine channels hg38 first and hg19 second
        bedpe_combined_ch = refor_hg38_ToBEDPE.out.filter{it[2]=='hg38'}.map{it->[it[1], it[0], it[2], it[3]]}.join(refor_hg19_ToBEDPE.out.filter{it[2]=='hg19'}.map{it->[it[1], it[0], it[2], it[3]]})

        //bedpe_combined_ch.view()
        
        emit:
            liftoverSVs_round2.out.filter{it[2]=='hg38'}.join(liftoverSVs.out)

    } else if (params.FromTo == 'hg38_to_hg19') {

        //run reformatToBEDPE on genomeTo

        reformatToBEDPE(crest_runs_ch.filter{it[2]=='hg19'})

        //pass only pedpe for hg38
        liftoverSVs(reformatToBEDPE.out.filter{it[2]=='hg38'})

        //combine channels hg38 first and hg19 second
        bedpe_combined_ch = liftoverSVs.out.map{it->[it[1], it[0], it[2], it[3]]}.join(reformatToBEDPE.out.filter{it[2]=='hg19'}.map{it->[it[1], it[0], it[2], it[3]]})

        emit:
            liftoverSVs.out.join(reformatToBEDPE.out.filter{it[2]=='hg19'}).view()
    } else {
        log.error("Not valid lift_over_fromTo. Valid modes are: 'GRCh37-lite|to|GRCh38' or 'GRCh38|to|GRCh37-lite' only ")
        exit 1
    }


    //pass pairs of bedpe (one liftover and one not) to compare the two using pairToPair  
    compareSVs(bedpe_combined_ch)

    //summarize 
    //summarizeOverlapResults(compareSVs.out)
    compareSVs.out.map{it->it[1]}
        .collectFile(name: "${params.outD}/bedpe_comp_${params.FromTo}/summary_all.out", keepHeader:true, skip:1, newLine:true)
        .view {file -> "pair\tNbreakpnts_hgFrom\tNbreakpnts_hgTo\tNoverlappingSVs\tNoverlap_match_type\tNoverlap_discordant_type\tNmissing\tNextra:${file.text}" }

    //get variant discordant rates for each pair for each sample used formula from functional equivalence paper (Nature Communincation 2018) 
    //formula for pair = (discordant + 0-only + 1-only + discordant_discordant_type)/(match+ discordant + match_discordant_type + discordant_discordant_type + 0-only + 1-only) 
    //match=overlapping positions AND same_strand AND sv-type is the same AND genotype is the same. match_discordant_type= all the same except sv-type. discordant=If there is overlap but genotypes are different irrespective if the sv-types are the same. discordant_discordant_type=if there is overlap and both genotype and sv-type are different.   
    
}//end of workflow

process compareSVs {
    label 'io_mem'
    publishDir "${params.outD}/bedpe_comp_${params.FromTo}", mode: 'copy'

    input:
    tuple val(pairName), val(anlsRunNameTo), val(genomeTo), file(bedpeFileTo), val(anlsRunNameFrom), val(genomeFrom), file(bedpeFileFrom)

    output:
    tuple val(pairName), file("summary_${pairName}_${params.FromTo}"), file("${pairName}_wliftover_both_samestrand.dbedpe"), file("${pairName}_wliftover_missing.bedpe"), file("${pairName}_wliftover_extra.bedpe"), emit: pairToPairOutput_ch

    """
#!/usr/bin/bash

module load bedtools/2.30.0
module load python/3.7.0

printf "pair\tNbreakpnts_hgFrom\tNbreakpnts_hgTo\tNoverlappingSVs\tNoverlapping_match_type\tNoverlapping_discordant_type\tEvidenceSummaryForOverlapping_diff\tEvidenceSummaryForOverlapping_mid\tNmissing\tEvidenceSummaryForMising\tNextra\tEvidenceSummaryForExtra\n" > summary_${pairName}_${params.FromTo}
echo -e "The From genome is the current genome and the To genome is the new genome \n"

##### only the first 10 columns matter
echo -e "The bedpeFileFrom can be two files depending on whether there was _inv file or not in the liftover step\n"
if [ `echo "${bedpeFileFrom}"|cut -f 1 -d " "|grep -v '_invs'|wc -l` == 1 ]; then
    targetBEDPE=`echo "${bedpeFileFrom}"|cut -f 1 -d " "`
else
    targetBEDPE=`echo "${bedpeFileFrom}"|cut -f 2 -d " "`
fi

echo -e "The From BEDPE file is set to \$targetBEDPE\n"

cat ${bedpeFileTo} |cut -f 1-10 > hgToF
echo -e "taking the From genome bedpe file with positions of breakpoint liftover to the To genome\n"
cat \$targetBEDPE |cut -f 1-10 > hgFromF
npnts_hgTo=`cat hgToF|wc -l`
npnts_hgFrom=`cat hgFromF|wc -l`

#### enforcing existence in both and ignoring strand

echo -e "option a is the From and option b is the To ==>finding overlapping features in both\n"
echo -e "will save overlapping features for pair to ${pairName}_wliftover_both_samestrand.dbedpe which is a double bedpe (20 columns) the first 10 are for the From and the next 10 are for the To\n"


echo -e "will use slop of ${params.slop} in the pairToPair to get overlaps" 
pairToPair  -type both -slop "${params.slop}" -is -a hgFromF -b hgToF | sort -k1,1V -k2,2n -k3,3n | uniq > ${pairName}_wliftover_both_samestrand.dbedpe

wait
##### added for testing, generally this gives the same results as the one above
pairToPair -type both -slop "${params.slop}"  -is -b hgFromF -a hgToF | sort -k1,1V -k2,2n -k3,3n | uniq > ${pairName}_wliftover_both_samestrand.dbedpe_flipped
######

wait

python ${workflow.projectDir}/bin/process_pairToPairOutput.py ${pairName}_wliftover_both_samestrand.dbedpe >Nmatch_discordant_types

nu_overlaps=`cat Nmatch_discordant_types| sed 1d|cut -f 1`
Noverlap_match_type=`cat Nmatch_discordant_types| sed 1d|cut -f 2`
Noverlap_discordant_type=`cat Nmatch_discordant_types| sed 1d|cut -f 3`
meanAbsDiffgeomCov=`cat Nmatch_discordant_types| sed 1d|cut -f 4`
meanMidGeomCov=`cat Nmatch_discordant_types| sed 1d|cut -f 5`

wait

## convert overlapping result to bedpe
cat ${pairName}_wliftover_both_samestrand.dbedpe | cut -f 1-10 > ${pairName}_wliftover_both_samestrand.bedpe_from
cat ${pairName}_wliftover_both_samestrand.dbedpe | cut -f 11-20 > ${pairName}_wliftover_both_samestrand.bedpe_to

wait

echo -e "finding SVs found using From genome but ends do not overlap the ends of previous overlapping result with the To genome (i.e. in From but exclude the overlap with the To)\n"
pairToPair  -type notboth -slop "${params.slop}" -is -a hgFromF -b ${pairName}_wliftover_both_samestrand.bedpe_from | sort -k1,1V -k2,2n -k3,3n | uniq > ${pairName}_wliftover_missing.bedpe

###### test if any of the missing are among the matched if there are then more work to be done to detect the missing
pairToPair  -type both -slop "${params.slop}" -is -a ${pairName}_wliftover_missing.bedpe -b ${pairName}_wliftover_both_samestrand.bedpe_from > test_anyMissing_stillInMatches

if [[ \$(cat test_anyMissing_stillInMatches|wc -l) -gt 0 ]]; then
    echo "still detecting some missing among the matches, they should be distinct ... exiting" > failures_detected
fi

wait
python ${workflow.projectDir}/bin/process_pairToPairOutput.py ${pairName}_wliftover_missing.bedpe >MeanCovACovB_missing
meangeomCovMissing=`cat MeanCovACovB_missing| sed 1d|cut -f 1`

echo -e "finding SVs found using To genome but ends do not overlap the ends of previous overlapping result (i.e. in To but exclude the overlap with the From)\n"
pairToPair -type notboth -slop "${params.slop}" -is -a hgToF -b ${pairName}_wliftover_both_samestrand.bedpe_to  | sort -k1,1V -k2,2n -k3,3n | uniq > ${pairName}_wliftover_extra.bedpe

###### test if any of the extra are among the matched if there are then more work to be done to detect the extra
pairToPair  -type both -slop "${params.slop}" -is -a ${pairName}_wliftover_extra.bedpe -b ${pairName}_wliftover_both_samestrand.bedpe_to > test_anyExtra_stillInMatches

if [[ \$(cat test_anyExtra_stillInMatches|wc -l) -gt 0 ]]; then
    echo "still detecting some extra among the matches, they should be distinct ... exiting" >> failures_detected
fi

wait
python ${workflow.projectDir}/bin/process_pairToPairOutput.py ${pairName}_wliftover_extra.bedpe >MeanCovACovB_extra
meangeomCovExtra=`cat MeanCovACovB_extra| sed 1d|cut -f 1`

wait

nu_missing=`cat ${pairName}_wliftover_missing.bedpe|wc -l`
nu_extra=`cat ${pairName}_wliftover_extra.bedpe|wc -l`

if [ ! -f failures_detected ]; then
    printf "${pairName}\\t\$npnts_hgFrom\\t\$npnts_hgTo\\t\$nu_overlaps\\t\$Noverlap_match_type\\t\$Noverlap_discordant_type\\t\$meanAbsDiffgeomCov\\t\$meanMidGeomCov\\t\$nu_missing\\t\$meangeomCovMissing\\t\$nu_extra\\t\$meangeomCovExtra" >>  summary_${pairName}_${params.FromTo}
else
    printf "${pairName}\\t\$npnts_hgFrom\\t\$npnts_hgTo\\t\$nu_overlaps\\t\$Noverlap_match_type\\t\$Noverlap_discordant_type\\t\$meanAbsDiffgeomCov\\t\$meanMidGeomCov\\t\failures_detected\\t\$meangeomCovMissing\\t\failures_detected\\t\$meangeomCovExtra" >>  summary_${pairName}_${params.FromTo}
fi

echo "done"

    """
    
}

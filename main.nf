#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl = 2

/* 
 * Run using the following command: 
 * 
 *   nextflow run reference_update_comparison_pipeline.nf 
 * 
 * 
 * The following parameters can be provided as command line options  
 * replacing the prefix `params.` with `--` e.g.:   
 * 
 *   nextflow run reference_update_comparison_pipeline.nf --sampleInfo ./sampleInfo.csv
 *   Ramzi Alsallaq
 * 
 */
/* set parameters to default values */
params.scratchDir = null 
params.mode = null // can be "crest" for crest output comparison, "crest-post" for crest-post comparison, or ... 

/* from official name of geneomeFrom to official name of genomeTo */
params.lift_over_fromTo='GRCh37-lite|to|GRCh38'

// crest params
params.crest_anls_run_names_file = null 
params.slop = 0 // wiggle in pairToPair

// crest-post params
params.crest_post_anls_run_names_file = null


params.outD = "$PWD/analysis/"

println("running mode = ",params.mode)

// default mode to use
def default_mode = "crest"

// variable to use throughout the workflow
def runMod = null

// check if CLI arg was passed; if so used that instead
if(params.mode == null){
    runMode = default_mode
    } else if(params.mode == "crest" | params.mode == "crest-post" | params.mode == "crest-localized" | params.mode == "crest-localized-post" | params.mode == "crest-tumoronly"| params.mode=="crest-tumoronly-post"| params.mode=="discordant" | params.mode=="discordant-post"| params.mode=="conserting-crest") {
        runMode = params.mode
} else {
    log.error("No valid mode specified. Valid modes are: 'crest', 'crest-post', ...etc ")
    exit 1
}

if (params.lift_over_fromTo == null) {
    log.error("not a valid value for lift_over_fromTo parameter, valid values are:GRCh37-lite|to|GRCh38 or GRCh38|to|GRCh37-lite")
    exit 1
} else if (params.lift_over_fromTo == 'GRCh37-lite|to|GRCh38'){
    params.FromTo = 'hg19_to_hg38'
} else if (params.lift_over_fromTo == 'GRCh38|to|GRCh37-lite') {
    params.FromTo = 'hg38_to_hg19'
} else {
    log.error("not a valid value for lift_over_fromTo parameter, valid values are:GRCh37-lite|to|GRCh38 or GRCh38|to|GRCh37-lite")
    exit 1
}


def helpMsg() {
log.info"""

Usage:
nextflow run main.nf <params>

reference_update_comparison_pipeline 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
runMode: $runMode
crest_anls_run_names_file : $params.crest_anls_run_names_file
crest_post_anls_run_names_file : $params.crest_post_anls_run_names_file
outD : $params.outD
scratchDir : $params.scratchDir
params.lift_over_fromTo : $params.lift_over_fromTo
FromTo : $params.FromTo
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
""".stripIndent()
}

params.help=false
if (params.help || params.keySet().size() == 0) {
    helpMsg()
    exit 0
}

// Make sure that --outD ends with "/"
if (!params.outD.endsWith("/")){
    outD_folder = params.outD.concat("/")
} else {
    outD_folder = params.outD
}


//import crest workflow
include { crest_wf } from './modules/crest_wf' 

//import crest_post workflow
include { crest_post_wf } from './modules/crest_post_wf' 

workflow {
    main:



    if(runMode == "crest") { // run the entire crest workflow 

        /* create channel from CSV file for crest runs*/
        crestRuns = channel.fromPath(params.crest_anls_run_names_file).splitCsv(header: true, sep: ",")
        //anls_run_name,pair_name,raw_somatic_sv_file
        crestRuns.map{ pair -> [ pair.anls_run_name,pair.pair_name, pair.genome, file(pair.raw_somatic_sv_file)] }.set{crest_runs_ch}

        crest_wf(crest_runs_ch)

    } else if (runMode == "crest-post") { //run the entire crest-post workflow

        /* create channel from CSV file for crest-post runs*/
        crestPostRuns = channel.fromPath(params.crest_post_anls_run_names_file).splitCsv(header: true, sep: ",")
        //officialGenome,anls_run_name,formal_name,pair_name,genome,fusionGenes_file
        crestPostRuns.map{ pair -> [ pair.officialGenome,pair.anls_run_name,pair.formal_name,pair.pair_name, pair.genome, file(pair.fusionGenes_file)] }.set{crest_post_runs_ch}

        crest_post_wf(crest_post_runs_ch)

   }
}

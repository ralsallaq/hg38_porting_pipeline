module load java/1.8.0_66
### get nextflow
#curl -s https://get.nextflow.io | bash
#####set up the environment
#source bin/set_env.sh

#### ./nextflow run reference_update_comparison_pipeline.nf --params.mode 'crest' --crest_anls_run_names_file ./anls_run_name.csv --outD ./analysis
 #./nextflow run main.nf --params.mode 'crest' --crest_anls_run_names_file ./run_test2.csv --outD ./analysis -resume
 ./nextflow run main.nf --mode 'crest-post-wcoord' --crest_post_anls_run_names_file "../anls_run_name.csv" --lift_over_fromTo 'GRCh37-lite|to|GRCh38' --slop 100 --outD ./analysis -resume

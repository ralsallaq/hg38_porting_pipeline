module load java/1.8.0_66
### get nextflow
#curl -s https://get.nextflow.io | bash
#####set up the environment
#source bin/set_env.sh

#### ./nextflow run reference_update_comparison_pipeline.nf --params.mode 'crest' --crest_anls_run_names_file ./anls_run_name.csv --outD ./analysis
 #./nextflow run main.nf --params.mode 'crest' --crest_anls_run_names_file ./run_test2.csv --outD ./analysis -resume
 #./nextflow run main.nf --mode 'crest-post-wcoord' --crest_post_anls_run_names_file "../anls_run_name.csv" --lift_over_fromTo 'GRCh37-lite|to|GRCh38' --slop 100 --outD ./analysis -resume

result_dir=/research_jude/rgs01_jude/groups/zhanggrp/projects/MethodDevelopment/common/ralsalla/minibam/optimized_branch/test_3_2023/frankenbam/recovery_events/test_manifest_frombed_bedpe/crest_hg38_ver_hg19/fastq_remapping_to_hg38
./nextflow run main.nf --params.mode 'crest'  --crest_anls_run_names_file $result_dir/crest_anls_run_names_file.csv --outD $result_dir/analysis_crest_comp -w $result_dir/work_crest_comp -resume

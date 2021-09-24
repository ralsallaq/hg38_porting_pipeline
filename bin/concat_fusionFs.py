import pandas as pd
manifest = pd.read_csv("anls_run_name.csv")
#sort so that files are aligned
manifest = manifest.sort_values(by='pair_name')
pairnames = manifest['pair_name'].drop_duplicates().values
hg19_dfs = [pd.read_csv(f, sep="\t") for f in manifest[manifest['genome']=='hg19']['fusionGenes_file'].values]
hg38_dfs = [pd.read_csv(f, sep="\t") for f in manifest[manifest['genome']=='hg38']['fusionGenes_file'].values]

hg19_all = pd.concat(hg19_dfs, axis=0)
hg38_all = pd.concat(hg38_dfs, axis=0)

hg19_all.to_csv("hg19_fusions_61samples.csv", index=False)
hg38_all.to_csv("hg38_fusions_61samples.csv", index=False)

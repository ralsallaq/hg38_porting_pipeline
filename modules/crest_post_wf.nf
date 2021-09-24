/*============================= crest-post workflow ===============================*/

include {concatFusions as concatFusions_hg19 } from './miscellaneous.nf'
include {concatFusions as concatFusions_hg38 } from './miscellaneous.nf'

workflow crest_post_wf {
    take:crest_post_runs_ch

    main:
    
    
    if ( params.FromTo == 'hg19_to_hg38') {

    concatFusions_hg19(crest_post_runs_ch.filter{it[4]=='hg19'}.map {it->it[5]}.toList(), channel.value('hg19'))
    concatFusions_hg38(crest_post_runs_ch.filter{it[4]=='hg38'}.map {it->it[5]}.toList(), channel.value('hg38'))

    predSVFusions_combined_ch = concatFusions_hg38.out.combine(concatFusions_hg19.out)

    //predSVFusions_combined_ch.view()
    compareFusionGenes(predSVFusions_combined_ch)

    
    }


}//end of workflow

process compareFusionGenes {
    label 'io_mem'
    publishDir "${params.outD}/crest-post/", mode: 'copy'

    input:
    tuple path(fusionGenesFileTo), path(fusionGenesFileFrom)

    output:
    path("summary_comp_fusionGenes.csv"), emit: fusionGeneCompOutput_ch

    """
#!/usr/bin/env python3
## define useful functions
def getNotNullFusionsForPair(pair):
    df19 = hg19_all[hg19_all['Sample(s)']==pair]
    df38 = hg38_all[hg38_all['Sample(s)']==pair]
    ### drop nulls
    df19 = df19[~df19['Fusion Gene'].isnull()]
    df38 = df38[~df38['Fusion Gene'].isnull()]
    cols = ['Sample(s)','ChrA','PosA','OrientA','GeneA','NumReadsA','ChrB', 'PosB', \
            'OrientB','GeneB','NumReadsB','Type','Usage','Fusion Gene','CDS','rating']
    df19 = df19[cols]
    df38 = df38[cols]
    return df38, df19

def meltGeneToSet(row):
    setofGenes = set(row.lower().split(","))
    return setofGenes

def findIndexOfMatch(row, fromGenome_df):
    for i, rr in fromGenome_df.iterrows():
        intsxSet = row['geneSet'].intersection(rr['geneSet'])
        if len(intsxSet)==0:
            continue
        else:
            return intsxSet, i
def compareFusionsForPair(pair, missingHQonly=False, extraHQonly=False):
    from itertools import product
    df_summary=pd.DataFrame()
    df38, df19 = getNotNullFusionsForPair(pair)
    Nfusions_hg38 = df38.shape[0] ## this number could include the same fusion genes occuring at different locations
    Nfusions_hg19 = df19.shape[0] ## this number could include the same fusion genes occuring at different locations
    
    ### sort for quicker access
    df38 = df38.sort_values(by="Fusion Gene")
    df19 = df19.sort_values(by="Fusion Gene")
    
    df38.loc[:,"geneSet"] = df38['Fusion Gene'].apply(lambda r: meltGeneToSet(r))
    df19.loc[:,"geneSet"] = df19['Fusion Gene'].apply(lambda r: meltGeneToSet(r))


    from itertools import product
    idx = 0
    for (i1, row1), (i2, row2) in product(df38.iterrows(),df19.iterrows()):
        #check for overlap between the two:
        intrsx = row1['geneSet'].intersection(row2['geneSet'])
        
        #skip rows that do not have overlaps
        if len(intrsx)==0:
            continue
        else:
#     cols = ['Sample(s)','ChrA','PosA','OrientA','GeneA','NumReadsA','ChrB', 'PosB', \
#             'OrientB','GeneB','NumReadsB','Type','Usage','Fusion Gene','CDS','rating']
            df_summary.loc[idx,'pair']=pair
            ### to be able to convert it to bedpe
            df_summary.loc[idx,'ChrA_hg38']=row1['ChrA']
            df_summary.loc[idx,'ChrA_hg19']=row2['ChrA']
            df_summary.loc[idx,'ChrB_hg38']=row1['ChrB']
            df_summary.loc[idx,'ChrB_hg19']=row2['ChrB']
            df_summary.loc[idx,'PosA_hg38']=row1['PosA']
            df_summary.loc[idx,'PosA_hg19']=row2['PosA']
            df_summary.loc[idx,'PosB_hg38']=row1['PosB']
            df_summary.loc[idx,'PosB_hg19']=row2['PosB']
            df_summary.loc[idx,'Type_hg38']=row1['Type']
            df_summary.loc[idx,'Type_hg19']=row2['Type']
            df_summary.loc[idx,'OrtA_hg38']=row1['OrientA']
            df_summary.loc[idx,'OrtA_hg19']=row2['OrientA']
            df_summary.loc[idx,'OrtB_hg38']=row1['OrientB']
            df_summary.loc[idx,'OrtB_hg19']=row2['OrientB']
            df_summary.loc[idx,'evidence_hg38']=row1['NumReadsA']+row1['NumReadsB']
            df_summary.loc[idx,'evidence_hg19']=row2['NumReadsA']+row2['NumReadsB']
            ######
            df_summary.loc[idx,'hg38Set'] = [row1['geneSet']]
            df_summary.loc[idx,'hg19Set'] = [row2['geneSet']]
            #print(intrsx)
            df_summary.loc[idx,'intrsx'] = [intrsx]
            df_summary.loc[idx,'hg38Rating'] = row1['rating']
            df_summary.loc[idx,'hg19Rating'] = row2['rating']
            df_summary.loc[idx,'hg38Usage'] = row1['Usage']
            df_summary.loc[idx,'hg19Usage'] = row2['Usage']
            idx +=1

    
    
    #### hold a set of genes from the overlap
    overlap_hold_set = set() 
    for i, row in df_summary.iterrows():
        if type(row["intrsx"]) is set:
            overlap_hold_set = overlap_hold_set.union(row["intrsx"])
        elif type(row["intrsx"]) is list:
            overlap_hold_set = overlap_hold_set.union(row["intrsx"][0])
    #print("the overlap set for ",pair," is ", len(overlap_hold_set))
            
    #### determine extra HQ
    for i1, row1 in df38.iterrows():
        ### in hg38 but not in the overlap
        extraSet = set(row1['geneSet'])-overlap_hold_set
        ### extra HQ
        if len(extraSet)>0 and row1['rating'] == "SJHQ" and extraHQonly:
            idx +=1
            df_summary.loc[idx,'extraHQ'] = [extraSet]
            df_summary.loc[idx,'pair']=pair
            df_summary.loc[idx,'hg38Set'] = [row1['geneSet']]
            df_summary.loc[idx,'hg38Rating'] = row1['rating']
            df_summary.loc[idx,'hg38Usage'] = row1['Usage']
            
            ### to be able to convert it to bedpe
            df_summary.loc[idx,'ChrA_hg38']=row1['ChrA']
            df_summary.loc[idx,'ChrB_hg38']=row1['ChrB']
            df_summary.loc[idx,'PosA_hg38']=row1['PosA']
            df_summary.loc[idx,'PosB_hg38']=row1['PosB']
            df_summary.loc[idx,'Type_hg38']=row1['Type']
            df_summary.loc[idx,'OrtA_hg38']=row1['OrientA']
            df_summary.loc[idx,'OrtB_hg38']=row1['OrientB']
            df_summary.loc[idx,'evidence_hg38']=row1['NumReadsA']+row1['NumReadsB']
            
            
        ### extra regardless to quality
        elif len(extraSet)>0 and not extraHQonly:
            idx +=1
            df_summary.loc[idx,'extra'] = [extraSet]
            df_summary.loc[idx,'pair']=pair
            df_summary.loc[idx,'hg38Set'] = [row1['geneSet']]
            df_summary.loc[idx,'hg38Rating'] = row1['rating']
            df_summary.loc[idx,'hg38Usage'] = row1['Usage']
            
            ### to be able to convert it to bedpe
            df_summary.loc[idx,'ChrA_hg38']=row1['ChrA']
            df_summary.loc[idx,'ChrB_hg38']=row1['ChrB']
            df_summary.loc[idx,'PosA_hg38']=row1['PosA']
            df_summary.loc[idx,'PosB_hg38']=row1['PosB']
            df_summary.loc[idx,'Type_hg38']=row1['Type']
            df_summary.loc[idx,'OrtA_hg38']=row1['OrientA']
            df_summary.loc[idx,'OrtB_hg38']=row1['OrientB']
            df_summary.loc[idx,'evidence_hg38']=row1['NumReadsA']+row1['NumReadsB']
            
        
    #### determine missing set
    for i2, row2 in df19.iterrows():
        ### in hg19 but not in the overlap
        missingSet = set(row2['geneSet'])-overlap_hold_set
        ### missing HQ
        if len(missingSet)>0 and row2['rating'] == "SJHQ" and missingHQonly:
            idx +=1
            df_summary.loc[idx,'missingHQ'] = [missingSet]
            df_summary.loc[idx,'pair']=pair
            df_summary.loc[idx,'hg19Set'] = [row2['geneSet']]
            df_summary.loc[idx,'hg19Rating'] = row2['rating']
            df_summary.loc[idx,'hg19Usage'] = row2['Usage']
            
            ### to be able to convert it to bedpe
            df_summary.loc[idx,'ChrA_hg19']=row2['ChrA']
            df_summary.loc[idx,'ChrB_hg19']=row2['ChrB']
            df_summary.loc[idx,'PosA_hg19']=row2['PosA']
            df_summary.loc[idx,'PosB_hg19']=row2['PosB']
            df_summary.loc[idx,'Type_hg19']=row2['Type']
            df_summary.loc[idx,'OrtA_hg19']=row2['OrientA']
            df_summary.loc[idx,'OrtB_hg19']=row2['OrientB']
            df_summary.loc[idx,'evidence_hg19']=row2['NumReadsA']+row2['NumReadsB']
            
        ### missing regardless to quality
        elif len(missingSet)>0 and not missingHQonly:
            idx +=1
            df_summary.loc[idx,'missing'] = [missingSet]
            df_summary.loc[idx,'pair']=pair
            df_summary.loc[idx,'hg19Set'] = [row2['geneSet']]
            df_summary.loc[idx,'hg19Rating'] = row2['rating']
            df_summary.loc[idx,'hg19Usage'] = row2['Usage']
            
            ### to be able to convert it to bedpe
            df_summary.loc[idx,'ChrA_hg19']=row2['ChrA']
            df_summary.loc[idx,'ChrB_hg19']=row2['ChrB']
            df_summary.loc[idx,'PosA_hg19']=row2['PosA']
            df_summary.loc[idx,'PosB_hg19']=row2['PosB']
            df_summary.loc[idx,'Type_hg19']=row2['Type']
            df_summary.loc[idx,'OrtA_hg19']=row2['OrientA']
            df_summary.loc[idx,'OrtB_hg19']=row2['OrientB']
            df_summary.loc[idx,'evidence_hg19']=row2['NumReadsA']+row2['NumReadsB']

            

    df_summary.loc[:,'hg19_totalFusions']  =  Nfusions_hg19
    df_summary.loc[:,'hg38_totalFusions']  =  Nfusions_hg38

    return df38, df19, df_summary

##### process the files

##tuple path(fusionGenesFileTo), path(fusionGenesFileFrom)
hg38_all = "${fusionGenesFileTo}"
hg19_all = "${fusionGenesFileFrom}"
print("heads of files")
print("The ToGenome file for combined fusions:\n")
print(hg38_all.head())
print("The FromGenome file for combined fusions:\n")
print(hg19_all.head())

print("tails of files")
print("The ToGenome file for combined fusions:\n")
print(hg38_all.tails())
print("The FromGenome file for combined fusions:\n")
print(hg19_all.tails())

print("processing the files\n")

pairsTo = hg38_all['Sample(s)'].sort_values().drop_duplicates().values
pairsFrom = hg19_all['Sample(s)'].sort_values().drop_duplicates().values

assert(np.all(pairsTo==pairsFrom)),"sample pairs are not aligned between the fusion results from the two genomes; check for missing data from some pairs in either genome"

df_summary_dfs = []
hg38_dfs = []
hg19_dfs = []
for p in pairsTo:
    print(p)
    temp1,temp2,temp3 = compareFusionsForPair(p)
    hg38_dfs.append(temp1)
    hg19_dfs.append(temp2)
    df_summary_dfs.append(temp3)

print("combine summaries from all pairs as is \n")
df_summary_all = pd.concat(df_summary_dfs, axis=0)

print("simplify how the genes are stored (strings)\n")
def convertRowToStr(row):
    if pd.isna(row):
        return None
    elif type(row) is set:
        return ",".join(row) #list(row)[0]
    elif type(row)is list:
        return ",".join(row[0])

def convertColToStr(col):
    df_summary_all.loc[:,col] = df_summary_all[col].apply(lambda r:convertRowToStr(r))

for c in ["hg38Set", "hg19Set","intrsx","missing","extra"]:
    convertColToStr(c)


print("save combined summaries into a CSV file\n")
df_summary_all.to_csv("summary_comp_fusionGenes.csv"

print("Done")

    """
}

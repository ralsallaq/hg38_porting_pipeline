# coding: utf-8
import pandas as pd
toF = pd.read_csv("hgToF", sep="\t", header=None)
fromF = pd.read_csv("hgFromF", sep="\t", header=None)
toF.head()
fromF.head()
def flipAdd(row):
    newRow = pd.Series(index=row.index.values)
    newRow.loc[0] = row.loc[3]
    newRow.loc[1] = row.loc[4]
    newRow.loc[2] = row.loc[5]
    newRow.loc[3] = row.loc[0]
    newRow.loc[4] = row.loc[1]
    newRow.loc[5] = row.loc[2]
    newRow.loc[6:7] = row.loc[6:7]
    newRow.loc[8] = '+' if row.loc[8] == '-' else '-'
    newRow.loc[9] = '+' if row.loc[9] == '-' else '-'
    
row = fromF.loc[0]
row
flipAdd(row)
def flipAdd(row):
    newRow = pd.Series(index=row.index.values)
    newRow.loc[0] = row.loc[3]
    newRow.loc[1] = row.loc[4]
    newRow.loc[2] = row.loc[5]
    newRow.loc[3] = row.loc[0]
    newRow.loc[4] = row.loc[1]
    newRow.loc[5] = row.loc[2]
    newRow.loc[6:7] = row.loc[6:7]
    newRow.loc[8] = '+' if row.loc[8] == '-' else '-'
    newRow.loc[9] = '+' if row.loc[9] == '-' else '-'
    return newRow 
    
flipAdd(row)
row
fromF_extended = fromF.copy()
additional=[]
for i, row in fromF.iterrows():
    additional.append(flipAdd(row))
    
additional
additional_df = pd.concat(additional, axis=0)
additional_df.head()
additional_df = pd.concat(additional, axis=1)
additional_df.head()
additional[0]
additional_df.head()
additional_df.T
fromF_extended.append(additional_df.T)
fromF.shape
fromF_extended.append(additional_df.T).to_csv("hgFromF_", index=False, header=None, sep="\t")
get_ipython().system('head hgFromF_')
get_ipython().system('wc -l hgFromF_')
toF_extended = toF.copy()
additional=[]
for i, row in toF.iterrows():
    additional.append(flipAdd(row))
    
additional_df = pd.concat(additional, axis=1)
additional_df.T.head()
toF_extended.append(additional_df.T).to_csv("hgToF_", index=False, header=None, sep="\t")
toF.shape
get_ipython().system('wc -l hgToF_')

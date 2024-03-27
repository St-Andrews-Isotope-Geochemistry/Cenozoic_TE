import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# mydifflib.py
from difflib import SequenceMatcher
from heapq import nlargest as _nlargest

def get_close_matches_indexes(word, possibilities, n=3, cutoff=0.6):
    """Use SequenceMatcher to return a list of the indexes of the best 
    "good enough" matches. word is a sequence for which close matches 
    are desired (typically a string).
    possibilities is a list of sequences against which to match word
    (typically a list of strings).
    Optional arg n (default 3) is the maximum number of close matches to
    return.  n must be > 0.
    Optional arg cutoff (default 0.6) is a float in [0, 1].  Possibilities
    that don't score at least that similar to word are ignored.
    """

    if not n >  0:
        raise ValueError("n must be > 0: %r" % (n,))
    if not 0.0 <= cutoff <= 1.0:
        raise ValueError("cutoff must be in [0.0, 1.0]: %r" % (cutoff,))
    result = []
    s = SequenceMatcher()
    s.set_seq2(word)
    for idx, x in enumerate(possibilities):
        s.set_seq1(x)
        if s.real_quick_ratio() >= cutoff and \
           s.quick_ratio() >= cutoff and \
           s.ratio() >= cutoff:
            result.append((s.ratio(), idx))

    # Move the best scorers to head of list
    result = _nlargest(n, result)

    # Strip scores for the best n matches
    return [x for score, x in result]


## import data
cwd=Path(os.getcwd())
data_path=cwd.parent.absolute()/"data"



foram_df=pd.read_csv(data_path/"foram_dataframe.csv", index_col=0)

foram_df_boron=foram_df.dropna(subset='B11', axis=0)


archive_df=pd.read_csv(data_path/"bigdf_240226.csv", index_col=0)
archive_df['time']=pd.to_datetime(archive_df['time'], 
                                        dayfirst=True)
archive_df=archive_df.reset_index(drop=True)


iso_to_iso_gas={'B11':'B11_No Gas', 
                'Mg24':'Mg24_No Gas', 
                'Li7':'Li7_No Gas', 
                'Na23':'Na23_No Gas', 
                'Al27':'Al27_No Gas', 
                'Mg25':'Mg25_No Gas', 
                'Ba138':'Ba138_No Gas', 
                'Mo92':'Mo92_No Gas', 
                'S32':'S32_48_O2', 
                'Co59': 'Co59_75_O2', 
                'U238':'U238_No Gas', 
                'Cd111':'Cd111_No Gas', 
                'Sr88': 'Sr88_No Gas'}





foram_agilent_link=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 'archive_idx', 'run_name', 'run_order', 'time', 'sample_name', 
                                        'rank', 'sum_sq_resid'])


for i, row in foram_df_boron.iterrows():
    sample_dict={}
    for key, val in iso_to_iso_gas.items():
        if pd.isnull(row[key]):
            continue
        sample_dict[val]=row[key]
    
    
    pivot_archive=archive_df.loc[archive_df['isotope_gas'].isin(list(sample_dict.keys()))]
    
    pivot_archive=pivot_archive.loc[~pivot_archive['sample_name'].str.contains('stg', case=False)]
    pivot_archive=pivot_archive.loc[~pivot_archive['sample_name'].str.contains('8301f', case=False)]
    pivot_archive.reset_index(inplace=True, names='archive_idx')
    pivot_archive=pd.pivot_table(pivot_archive, columns='isotope_gas', 
                                 values='cali_single', 
                                 index=['run_name', 'run_order', 'sample_name', 'time', 'archive_idx'])
    pivot_archive.reset_index(inplace=True)
    
    sample_name_idx=np.array([False]*pivot_archive.shape[0], dtype=bool)
    sample_name_idx_2=np.array([False]*pivot_archive.shape[0], dtype=bool)
    
    #first search for sample_id and sample_mix_id in the archive sample_name
    if  pd.notnull(row['sample_id']):
        sample_name_idx=pivot_archive['sample_name'].str.contains(row['sample_id'], case=False)       
        if sum(sample_name_idx)>1:
            sample_name_idx=np.array([False]*pivot_archive.shape[0], dtype=bool)
            matches=get_close_matches_indexes(row['sample_id'], pivot_archive['sample_name'].values, n=1)
            sample_name_idx[matches[0]]=True

    if  pd.notnull(row['sample_mix_id']):
        sample_name_idx_2=pivot_archive['sample_name'].str.contains(row['sample_mix_id'], case=False)
        if sum(sample_name_idx_2)>1:
            sample_name_idx_2=np.array([False]*pivot_archive.shape[0], dtype=bool)
            matches=get_close_matches_indexes(row['sample_mix_id'], pivot_archive['sample_name'].values, n=1)
            sample_name_idx_2[matches[0]]=True 

    
    sample_name_idx=sample_name_idx | sample_name_idx_2
    s_df=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 
                               'archive_idx', 'run_name', 'run_order', 'time', 'sample_name', 
                               'rank', 'sum_sq_resid'])
    
    if sample_name_idx.sum()==1:
        s_df['archive_idx']=pivot_archive.loc[sample_name_idx, 'archive_idx']
        s_df['foram_database_idx']=i
        s_df['core']=row['core']
        s_df['sample_id']=row['sample_id']
        s_df['sample_mix_id']=row['sample_mix_id']
        s_df['run_name']=pivot_archive.loc[sample_name_idx, 'run_name']
        s_df['run_order']=pivot_archive.loc[sample_name_idx, 'run_order']
        s_df['time']=pivot_archive.loc[sample_name_idx, 'time']
        s_df['sample_name']=pivot_archive.loc[sample_name_idx, 'sample_name']
        s_df['rank']=0
        s_df['sum_sq_resid']=0
        s_df.set_index('foram_database_idx', inplace=True)
        foram_agilent_link=pd.concat([foram_agilent_link, s_df], axis=0)
        continue
    
    resid_dict={}
    
    #find the squared relative residuals for each element
    for key, val in sample_dict.items():

        sq_resid=(1-pivot_archive[key].values/val)**2
        
        resid_dict[key]=sq_resid
    
    resid_df=pd.DataFrame(resid_dict)
    
    resid_df['sum_sq_resid']=resid_df.sum(axis=1)
    
    resid_df.sort_values('sum_sq_resid', inplace=True)
    
    idx_list=resid_df.index[:10]
    
    s_df[['run_name', 'run_order', 'time', 'sample_name']]=pivot_archive.loc[idx_list, ['run_name', 'run_order', 'time', 'sample_name']].values
    s_df[['rank', 'archive_idx', 'sum_sq_resid']]=np.array([np.arange(1, 11), pivot_archive.loc[idx_list, 'archive_idx'].values, 
                                                            resid_df.loc[idx_list, 'sum_sq_resid'].values]).T
    s_df[['sample_id', 'sample_mix_id', 'core', 'foram_database_idx']]=[row['sample_id'], row['sample_mix_id'], row['core'], i]
    s_df.set_index('foram_database_idx', inplace=True)
    foram_agilent_link=pd.concat([foram_agilent_link, s_df], axis=0)


foram_agilent_link.to_csv(data_path/"foram_agilent_link.csv")
    

        
        





import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import duckdb

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
data_path=Path(os.getcwd())/"data"




foram_df=pd.read_csv(data_path/"foram_dataframe_240625.csv", index_col=0)

foram_df_boron=foram_df.dropna(subset='B11', axis=0)

foram_df_boron['s_id_species']=foram_df_boron['sample_mix_id']+'_'+foram_df_boron['species']




archive_df=pd.read_csv(data_path/"bigdf_240918.csv", index_col=0)
archive_df['time']=pd.to_datetime(archive_df['time'], 
                                      dayfirst=True)


archive_df=archive_df.reset_index(drop=True)


archive_df['agilent_id']='STG'+archive_df['time'].dt.strftime('%Y%m%d%H%M')

#remove duplicates
archive_df=archive_df.drop_duplicates(subset=['agilent_id', 'isotope_gas', 'sample_name', 'time', 'run_name', 'run_order'], keep='first')

archive_df.to_csv(data_path/"archive_cleaned_240923.csv")

# Connect to an in-memory DuckDB database
conn = duckdb.connect(database=':memory:')

#filter runs older than 2021
query = f"""
SELECT *
FROM read_csv_auto('{data_path/"bigdf_240918.csv"}')
WHERE time < '2021-01-01' AND time > '2019-12-12'
"""

# Query the CSV file directly without loading it into memory
query = f"""
SELECT *
FROM read_csv_auto('{data_path/"bigdf_240226.csv"}')
"""


query = f"""
SELECT DISTINCT sample_name
FROM read_csv_auto('{data_path/"bigdf_240226.csv"}')
WHERE run_name = 'Carbonates_TE_20smps_JB_20211215_2.b'
"""
# Execute the query and fetch the results into a pandas DataFrame
result = conn.execute(query).fetchdf()

core_runame_dict={
    '522':['Carbonates_TE_20smps_JB_20211215_2.b', 
           'Carbonates_TE_20smps_JB_20220218.b', 
           'Carbonates_TE_19smps_JB_RG_20230717.b', 
           'Carbonates_TE_22smps_JB_RG_20230803_0.5mmol.b', 
           'Carbonates_TE_22smps_JB_RG_20230807_2.b'], 
    '689B':['Carbonates_TE_18smps_JB_JSE_20230613_2.b', 
           'Carbonates_TE_20smps_JB_JSE_20230622.b', 
           'Carbonates_TE_22smps_JB_JSE_20230704.b', 
           'Carbonates_TE_22smps_JB_RG_20230807_2.b', 
           'Carbonates_TE_19smps_JB_GW_20230811.b'], 
    '1209':['Carbonates_TE_20smps_JB_20230420.b', 
            'Carbonates_TE_24smps_JB_20230426.b', 
            'Carbonates_TE_21smps_JB_20230504.b', 
            'Carbonates_TE_14smps_JB_20230718.b'], 
    'U1409':['Carbonates_TE_20smps_JB_GW_20230616.b', 
             'Carbonates_TE_21smps_JB_GW_20230808.b', 
             'Carbonates_TE_21smps_JB_GW_20230808.b', 
             'Carbonates_TE_19smps_JB_GW_20230811.b'], 
    '1264':['Carbonates_TE_21smps_JB_ED_20221026.b', 
            'Carbonates_TE_24smps_JB_ED_20221103.b', 
            'Carbonates_TE_13smps_HB_ED_20221117_V2.b',
            'Carbonates_TE_24smps_JB_ED_20221124.b', 
            'Carbonates_TE_14smps_JB_ED_20221212_2.b', 
            'Carbonates_TE_24smps_HB_ED_20221206.b'], 
    '1262':['Carbonates_TE_3smps_JB_20220428_3.b', 
            'Carbonates_TE_4smps_JB_20210408.b', 
            'Carbonates_TE_7smps_JB_20220428_4_1.b', 
            'Carbonates_TE_8smps_JB_20220428.b', 
            'Carbonates_TE_8smps_JB_20220428_2.b', 
            'Carbonates_TE_12smps_JB_20220302.b', 
            'Carbonates_TE_16smps_JB_20200203_#2.b', 
            'Carbonates_TE_17smps_JB_MS_20210408.b', 
            'Carbonates_TE_17smps_JB_20220914.b', 
            'Carbonates_TE_48smps_JB_MS_20210512.b', 
            'Carbonates_TE_14smps_JB_MS_20210617.b', 
            'Carbonates_TE_16smps_JB_20200203_#2.b', 
            'Carbonates_TE_22smps_JB_20211206.b', 
            'Carbonates_TE_54smps_JB_20210409.b', 
            'Carbonates_TE_28smps_JB_20210331_new.b', 
            'Carbonates_TE_28smps_JB_20210406.b']
}

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



foram_agilent_link=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 's_id_species', 'agilent_id', 'run_name', 'run_order', 'time', 'sample_name', 
                                        'rank', 'sum_sq_resid'])


# first check for exact matches by name


for i, row in foram_df_boron.iterrows():
    sample_dict={}
    for key, val in iso_to_iso_gas.items():
        if pd.isnull(row[key]):
            continue
        sample_dict[val]=row[key]
    
    core=row['core']
    
    pivot_archive=archive_df.loc[archive_df['isotope_gas'].isin(list(sample_dict.keys()))]
    pivot_archive=pivot_archive.loc[~pivot_archive['sample_name'].str.contains('stg', case=False)]
    pivot_archive=pivot_archive.loc[~pivot_archive['sample_name'].str.contains('8301f', case=False)]
    #pivot_archive.reset_index(inplace=True, names='archive_idx')
    pivot_archive=pd.pivot_table(pivot_archive, columns='isotope_gas', 
                                 values='cali_single', 
                                 index=['run_name', 'run_order', 'sample_name', 'time', 'agilent_id'])
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
    s_df=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 's_id_species', 
                               'agilent_id', 'run_name', 'run_order', 'time', 'sample_name', 
                               'rank', 'sum_sq_resid'])
    
    if sample_name_idx.sum()==1:
        s_df['agilent_id']=pivot_archive.loc[sample_name_idx, 'agilent_id']
        s_df['foram_database_idx']=i
        s_df['core']=row['core']
        s_df['sample_id']=row['sample_id']
        s_df['sample_mix_id']=row['sample_mix_id']
        s_df['s_id_species']=row['s_id_species']
        s_df['run_name']=pivot_archive.loc[sample_name_idx, 'run_name']
        s_df['run_order']=pivot_archive.loc[sample_name_idx, 'run_order']
        s_df['time']=pivot_archive.loc[sample_name_idx, 'time']
        s_df['sample_name']=pivot_archive.loc[sample_name_idx, 'sample_name']
        s_df['rank']=0
        s_df['sum_sq_resid']=0
        s_df.set_index('foram_database_idx', inplace=True)
        foram_agilent_link=pd.concat([foram_agilent_link, s_df], axis=0)
        

foram_agilent_link.drop_duplicates(inplace=True)



foram_agilent_link.to_csv(data_path/"exact_agilent_matches_240923.csv")
    
exact_matches=foram_agilent_link.copy()

isotope=['B11', 'Mg24', 'Al27', 'Na23', 'Li7', 'Mg25', 'Ba138',  'U238', 'Cd111', 'Sr88', 'Nd146', 'Mn55']


for iso in isotope:
    exact_matches[iso+'_foram']=foram_df_boron.loc[exact_matches.index, iso]




for i, row in exact_matches.iterrows():
    sample=archive_df.loc[archive_df['agilent_id']==row['agilent_id']]
    
    for iso in isotope:
        val=sample.loc[sample['isotope_gas']==iso+'_No Gas', 'cali_single'].values
        exact_matches.loc[i, iso+'_agilent']=val


import statsmodels.api as sm 
fit_dict={}
for iso in isotope: 
    df=exact_matches.dropna(subset=[iso+'_foram', iso+'_agilent'], axis=0)   
    x=sm.add_constant(df[iso+'_foram'])
    y=df[iso+'_agilent']
    fit=sm.OLS(y, x).fit()
    fit_dict[iso+'_No Gas']=fit
    print(iso)
    print(fit.rsquared)







#filter out the foram archive to remove exact matches
foram_df_unmatched=foram_df_boron.loc[~foram_df_boron.index.isin(exact_matches.index)]
#filter out the archive to remove exact matches
archive_df_unmatched=archive_df.loc[~archive_df['agilent_id'].isin(exact_matches['agilent_id'])]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('stg', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('8301f', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('stgcs', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('blk', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('cs1', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('8301c', case=False)]
archive_df_unmatched_B=archive_df_unmatched.loc[archive_df_unmatched['isotope_gas'].str.contains('B11_No Gas')]

archive_df_unmatched_B[['run_name', 'run_order', 'sample_name', 'time', 'agilent_id', 'brkt_stnd', 'isotope_gas', 'cali_single']].to_csv(data_path/"archive_df_unmatched_B.csv", index=False)

foram_df_unmatched.to_csv(data_path/"foram_dataframe_240923_unmatched.csv")



#cycle through cores and match the samples
link_dict={}


for key, val in core_runame_dict.items():
    print(key)
    foram_core_df=foram_df_unmatched.loc[foram_df_unmatched['core']==key]
    
    archive_df_core=archive_df_unmatched.loc[archive_df_unmatched['run_name'].isin(val)]
    
    if key == '1262':
        old_archive=archive_df_core.loc[archive_df_core['time'] < '2021-01-01']
        archive_df_core=pd.concat([archive_df_core, old_archive], axis=0)
    

    foram_agilent_link=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 's_id_species', 'agilent_id', 'run_name', 'run_order', 'time', 'sample_name', 
                                        'rank', 'sum_sq_resid'])

    
    for i, row in foram_core_df.iterrows():
        sample_dict={}
        for k, v in iso_to_iso_gas.items():
            if pd.isnull(row[k]):
                continue
            if v in fit_dict.keys():
    
                iso_vals=row[k]
                new_iso_vals=fit_dict[v].predict([1, iso_vals])
                
                sample_dict[v]=new_iso_vals
        
        
        pivot_archive=archive_df_core.loc[archive_df_core['isotope_gas'].isin(list(sample_dict.keys()))]
        #filter archive
        sample_mix_id=row['sample_mix_id']
        if '2020' in sample_mix_id:
            pivot_archive=pivot_archive.loc[pivot_archive['time'] < '2021-01-01']
        
        
        pivot_archive=pd.pivot_table(pivot_archive, columns='isotope_gas', 
                                    values='cali_single', 
                                    index=['run_name', 'run_order', 'sample_name', 'time', 'agilent_id'])

        
        pivot_archive.reset_index(inplace=True)
        
        
        sample_name_idx=np.array([False]*pivot_archive.shape[0], dtype=bool)
        sample_name_idx_2=np.array([False]*pivot_archive.shape[0], dtype=bool)
        
        
        sample_name_idx=sample_name_idx | sample_name_idx_2
        s_df=pd.DataFrame(columns=['foram_database_idx', 'core','sample_id', 'sample_mix_id', 's_id_species', 
                               'agilent_id', 'run_name', 'run_order', 'time', 'sample_name', 
                               'rank', 'sum_sq_resid'])
        
        resid_dict={}
        
        #find the squared relative residuals for each element
        for k, v in sample_dict.items():

            sq_resid=(1-pivot_archive[k].values/v)**2
            
            resid_dict[k]=sq_resid
        
        resid_df=pd.DataFrame(resid_dict)
        
        resid_df['sum_sq_resid']=resid_df.sum(axis=1)
        
        resid_df.sort_values('sum_sq_resid', inplace=True)
        
        idx_list=resid_df.index[:10]
        
        s_df[['run_name', 'run_order', 'time', 'sample_name']]=pivot_archive.loc[idx_list, ['run_name', 'run_order', 'time', 'sample_name']].values
        s_df[['rank', 'agilent_id', 'sum_sq_resid']]=np.array([np.arange(1, 11), pivot_archive.loc[idx_list, 'agilent_id'].values, 
                                                                resid_df.loc[idx_list, 'sum_sq_resid'].values]).T
        s_df[['sample_id', 'sample_mix_id', 's_id_species', 'core', 'foram_database_idx']]=[row['sample_id'], row['sample_mix_id'], row['s_id_species'], row['core'], i]
        s_df.set_index('foram_database_idx', inplace=True)
        foram_agilent_link=pd.concat([foram_agilent_link, s_df], axis=0)
        
        
    
    link_dict[key]=foram_agilent_link.copy()    
    foram_agilent_link.to_csv(data_path/f"link_{key}.csv")



exact_matches['match_type']='exact_name'
exact_matches['sum_sq_resid']=np.nan


matched_df=exact_matches.loc[:, ['sample_id', 'sample_mix_id', 's_id_species', 'agilent_id', 'match_type', 'sum_sq_resid']]


bad_matches=pd.DataFrame()
for key, val in link_dict.items():
    print(key)

    core_df=val.drop('foram_database_idx', axis=1).reset_index(names='foram_database_idx')
    core_matches=core_df.loc[core_df['rank']==1]

    duplicates=len(core_matches)-len(pd.unique(core_matches['agilent_id']))
    
    while duplicates>0:
        print(duplicates)    
        for unique_sample in pd.unique(core_matches['agilent_id']):

            
            #how many matches
            match_number=len(core_matches.loc[core_matches['agilent_id']==unique_sample])
            
            
            
            while match_number>1:

                idx=core_matches.loc[core_matches['agilent_id']==unique_sample, 'sum_sq_resid'].idxmax()

                rank=core_matches.loc[idx, 'rank']

                
                
                s_df=core_df.loc[core_df['s_id_species']==core_matches.loc[idx, 's_id_species']]
                
                if rank == 10:
                    bad_matches = pd.concat([bad_matches, s_df.loc[s_df['rank']==10]], axis=0)
                    core_matches.drop(idx, axis=0, inplace=True)
                else:
                    s=s_df.loc[s_df['rank']==rank+1]
                    
                    core_matches.loc[idx]=s.iloc[0, :]
                    
                match_number=len(core_matches.loc[core_matches['agilent_id']==unique_sample])

        duplicates=len(core_matches)-len(pd.unique(core_matches['agilent_id']))
        
    core_matches['match_type']='closest_match'
    matched_df=pd.concat([matched_df, core_matches.loc[:, ['sample_id', 'sample_mix_id', 's_id_species', 'agilent_id', 'match_type', 'sum_sq_resid']]], axis=0)


bad_matches['match_type']='over_rank'

bad_matches=bad_matches.loc[:, ['sample_id', 'sample_mix_id', 's_id_species', 'agilent_id', 'match_type', 'sum_sq_resid']]

bad_matches=pd.concat([bad_matches, matched_df.loc[matched_df['sum_sq_resid']>2]], axis=0)  
good_matches=matched_df.loc[(matched_df['sum_sq_resid']<=2)| (matched_df['match_type']=='exact_name')]




foram_df_boron_matched=foram_df_boron.copy()

isotopes=['B11', 'Mg24', 'Al27', 'Na23', 'Li7', 'Mg25', 'Ba138',  'U238', 'Cd111', 'Sr88', 'Nd146', 'Mn55']

for iso in isotopes:
    iso_archive=archive_df.loc[archive_df['isotope_gas'].str.contains(iso)]

    for i, row in good_matches.iterrows():
        foram_df_idx=foram_df_boron_matched['s_id_species'] == row['s_id_species']
        
        foram_df_boron_matched.loc[foram_df_idx, iso]=iso_archive.loc[iso_archive['agilent_id']==row['agilent_id'], 'cali_single'].values[0]
        
        



foram_df_boron_matched.to_csv(data_path/"foram_dataframe_240625_matched.csv")       
    
    

foram_matched_subset=foram_df_boron_matched[['s_id_species', 'species', 'species_simple']+isotopes].copy()

#get counts of nan values in each column
foram_matched_subset.isnull().sum()
    
foram_matched_subset.dropna(subset=['Mg24', 'Li7', 'Mn55'], axis=0, inplace=True)  
foram_matched_subset.drop(['Mg25', 'Nd146'], axis=1, inplace=True)    
    
foram_matched_subset_noU=foram_matched_subset.dropna(subset=['U238'], axis=0)   



foram_matched_subset_noU.to_csv(data_path/"foram_dataframe_240625_matched_noU.csv")








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




archive_df=pd.read_csv(data_path/"archive_cleaned_240923.csv", index_col=0)
archive_df['time']=pd.to_datetime(archive_df['time'], 
                                      dayfirst=True)
archive_df=archive_df.reset_index(drop=True)





# Connect to an in-memory DuckDB database
conn = duckdb.connect(database=':memory:')

#filter runs older than 2021
query = f"""
SELECT *
FROM read_csv_auto('{data_path/"bigdf_240226.csv"}')
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




#get up to date matches
foram_agilent_link=pd.read_csv(data_path/"exact_agilent_matches.csv", index_col=0)
James_matching_df=pd.read_csv(data_path/"foram_dataframe_240625_unmatched_JB.csv", index_col=0)


JB_questions=James_matching_df.dropna(subset=['matching_comment'], axis=0)
JB_matches=James_matching_df.dropna(subset=['agilent_id'], axis=0)


foram_df_w_agilent_id=pd.merge(foram_df, foram_agilent_link['agilent_id'], how='left', 
                               left_index=True, right_index=True)



foram_df_w_agilent_id.loc[JB_matches.index, 'agilent_id']=JB_matches['agilent_id']



foram_df_w_agilent_id_boron=foram_df_w_agilent_id.dropna(subset='B11', axis=0)


unmatched_foram_B=foram_df_w_agilent_id_boron[pd.isnull(foram_df_w_agilent_id_boron['agilent_id'])]





## matching by inference




sources_to_ignore=['Greenop']


unmatched_foram_B_source_filtered=unmatched_foram_B[~unmatched_foram_B['source'].isin(sources_to_ignore)]

archive_df_unmatched=archive_df[~archive_df['agilent_id'].isin(foram_df_w_agilent_id_boron['agilent_id'])]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('stg', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('8301f', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('stgcs', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('blk', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('cs1', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('8301c', case=False)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['sample_name'].str.contains('wash', case=False)]


archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['run_name'].str.contains('STDS', case=True)]
archive_df_unmatched=archive_df_unmatched.loc[~archive_df_unmatched['run_name'].str.contains('_MD_', case=True)]



Elena_runs=archive_df_unmatched.loc[archive_df_unmatched['run_name'].str.contains('ED', case=True), 'run_name'].unique()
Wayne_runs=pd.unique(archive_df_unmatched.loc[archive_df_unmatched['run_name'].str.contains('WS', case=True), 'run_name'])

archive_df_Elena=archive_df_unmatched.loc[archive_df_unmatched['run_name'].isin(Elena_runs)]
archive_df_Wayne=archive_df_unmatched.loc[archive_df_unmatched['run_name'].isin(Wayne_runs)]


archive_df_Elena_B=archive_df_Elena.loc[archive_df_Elena['isotope_gas']=='B11_No Gas', 'cali_single'].values
archive_df_Elena_Mg=archive_df_Elena.loc[archive_df_Elena['isotope_gas']=='Mg24_No Gas', 'cali_single'].values
unmatched_foram_B_Elena=unmatched_foram_B_source_filtered[unmatched_foram_B_source_filtered['source']=='Elena']


fig, ax =plt.subplots()
plt.scatter(archive_df_Elena_B, archive_df_Elena_Mg, label='archive')
plt.scatter(unmatched_foram_B_Elena['B11'], unmatched_foram_B_Elena['Mg24'], label='unmatched', s=20)
plt.legend()









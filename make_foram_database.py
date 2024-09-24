
from pathlib import Path
import os
import pandas as pd
import numpy as np
import time

def depth_m_to_pressure_bar(depth):
    pressure_bar=1023.6*9.80665*depth*10**-5
    return pressure_bar

#paths
data_path=Path(os.getcwd())/"data"
fig_path=Path(os.getcwd())/"figures"




## Foram data   
#import foram data
excel_file_path = data_path/"foram_database_unstructured.xlsx"

# Read all sheets from the Excel file into a dictionary of DataFrames
excel_data = pd.read_excel(excel_file_path, sheet_name=None)

foram_df=pd.DataFrame()
# Iterate through the sheets and store each DataFrame in the dictionary
for sheet_name, core_df in excel_data.items():
    core_df.insert(0, 'core', sheet_name)
    foram_df=pd.concat([foram_df, core_df])

foram_df.reset_index(drop=True, inplace=True)


#set lat lons and basin
lats={'1262':-27.186, '522':-26.114, '689B': -64.517000, '1264': -28.532680, '999':12.744, 
      '926':3.719017, '1209':32.6517, 'U1409_U1407':41.42498833, '1218':8.889630, 'U1334':7.999, 
      '1263':-28.532830, 'U1406': 40.349930}
lons={'1262':1.577, '522':-5.130, '689B':3.0999, '1264': 2.845510, '999': -78.7393, 
      '926':-42.908300, '1209':158.505983, 'U1409_U1407':-49.8133117, '1218': -135.366660, 'U1334':-131.973, 
      '1263':2.779480, 'U1406': -51.649830}
basin={'1262':'Atlantic', '522':'Atlantic', '689B':'Southern', '1264':'Atlantic', '999': 'Atlantic', 
       '926':'Atlantic', '1209':'Pacific', 'U1409_U1407':'Atlantic', '1218':'Pacific', 'U1334':'Pacific', 
       '1263':'Atlantic', 'U1406': 'Atlantic'}

foram_df['lat']=foram_df['core'].map(lats)
foram_df['lon']=foram_df['core'].map(lons)
foram_df['ocean_basin']=foram_df['core'].map(basin)

#1262 paleodepth
foram_df.loc[foram_df['core']=='1262', 'palaeo_depth_m']=np.interp(foram_df.loc[foram_df['core']=='1262', 'age_Ma'], [56, 66], [3000, 3500])


#convert sample names that are single digits to core_digit.
foram_df['sample_id']=np.where(foram_df['sample_id'].str.isdigit(), 
                               foram_df['core'].astype(str)+'_'+foram_df['sample_id'].astype(str), 
                               foram_df['sample_id'])

foram_df['sample_mix_id']=np.where(foram_df['sample_id'].str.isdigit(), 
                               foram_df['core'].astype(str)+'_'+foram_df['sample_mix_id'].astype(str), 
                               foram_df['sample_mix_id'])

#make sure species names are consistent
foram_df.loc[foram_df['species']== 'Cibicidoides (prae)mundulus', 'species']= 'Cibicidoides praemundulus'


#make a simplified species column to account for different mundulus names
foram_df['species_simple']=foram_df['species']
foram_df.loc[foram_df['species'].str.contains('mundulus'), 'species_simple']='Cibicidoides mundulus'
foram_df.loc[foram_df['species'].str.contains('eocaenus'), 'species_simple']='Cibicidoides eocaenus'
foram_df.loc[foram_df['species'].str.contains('grimsdalei'), 'species_simple']='Cibicidoides grimsdalei'



#convert depth to pressure
foram_df['pressure_bar']=depth_m_to_pressure_bar(foram_df['palaeo_depth_m'].values)


timestr = time.strftime("%y%m%d-%H%M")
#save to disk
foram_df.to_csv(data_path/f"foram_database_{timestr}.csv")

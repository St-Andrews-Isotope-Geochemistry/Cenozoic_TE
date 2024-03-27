from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


data_path=Path(os.getcwd())/"data"

#import data
foram_df=pd.read_csv(data_path/"foram_dataframe.csv", index_col=0)


#subset data to only include TE data (using B11 as a stand in for TE data)
foram_te_df=foram_df.loc[pd.notnull(foram_df['B11'])]


#what is the distribution of the species through time?
fig, ax =plt.subplots(figsize=(10, 6))
sns.stripplot(data=foram_te_df, x='age_Ma', hue='species_simple', 
              size=15, y='species_simple', ax=ax, alpha=.1, 
              jitter=False, linewidth=1, marker="s")
ax.set(ylabel=None)



#Given the distributions, let's omit velascoensis, hyphalus, umbonifera, beccariiformis, subspiratus, and laurisae
#omit_species=['velascoensis', 'hyphalus', 'umbonifera', 'beccariiformis', 'laurisae', 'subspiratus']
#for species in omit_species:
#    foram_te_df=foram_te_df.loc[~foram_te_df['species_simple'].str.contains(species)]

isotopes=['Li7', 'B11', 'Na23', 'Mg24', 'Al27', 'P31', 'S32', 'K39',
      'Mn55','Fe56','Ni58','Co59','Cu63', 'Zn64', 'Rb85', 'Sr88', 
      'Mo92','Cd111', 'Ba138', 'Nd146', 'U238']


for iso in isotopes:
    fig, ax =plt.subplots(figsize=(12, 6))
    sns.scatterplot(data=foram_te_df, x='age_Ma', y=iso, hue='species_simple')
    plt.legend(fontsize=6, title=None)
    plt.ylabel(iso+'/Ca')


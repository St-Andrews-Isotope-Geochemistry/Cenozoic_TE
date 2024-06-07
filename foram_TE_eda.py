from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



def find_outliers(array1, mod=1.5):

    """
    Returns boolean array where true values denote outliers in original array
    
    Parameters
    ----------
    array1 : 1d or 2d array (numpy array)
        array to search for outliers.
    mod : modifier of outlier distance (iqr multiplier), default 1.5.

    Returns
    -------
    retarray : numpy array of len(array1)
        An array of boolean values where True values indicate the presence of 
        outliers at the same index of array1.

    """
    array1=np.array(array1, dtype=float)
    array1=array1.flatten()
    x = array1[~np.isnan(array1)]
    if len(x)>2:
        q75, q25 = np.percentile(x, [75 ,25])
        iqr = q75 - q25
        outs=((array1>iqr*mod+q75) | (array1<q25-iqr*mod))
    else:
        outs=np.isnan(array1)
        
    return outs
def isotope_to_element(isotope):
    """
    Returns the element of the isotope
    
    Parameters
    ----------
    isotope : str
        The name of the isotope.

    Returns
    -------
    str
        The element of the isotope.

    """
    elements={'Li7':'Li', 'B11':'B', 'Na23':'Na', 'Mg24':'Mg', 'Al27':'Al', 'P31':'P', 'S32':'S', 'K39':'K',
      'Mn55':'Mn','Fe56':'Fe','Ni58':'Ni','Co59':'Co','Cu63':'Cu', 'Zn64':'Zn', 'Rb85':'Rb', 'Sr88':'Sr', 
      'Mo92':'Mo','Cd111':'Cd', 'Ba138':'Ba', 'Nd146':'Nd', 'U238':'U'}
    return elements[isotope]
def isotope_to_units(isotope):
    """
    Returns the units of the isotope
    
    Parameters
    ----------
    isotope : str
        The name of the isotope.

    Returns
    -------
    str
        The units of the isotope.

    """
    units={'Li7':r'$\mu$mol/mol', 'B11':r'$\mu$mol/mol', 'Na23':'mmol/mol', 'Mg24':'mmol/mol', 
           'Al27':r'$\mu$mol/mol', 'P31':r'$\mu$mol/mol', 'S32':'mmol/mol', 'K39':r'$\mu$mol/mol',
            'Mn55':r'$\mu$mol/mol', 'Fe56':r'$\mu$mol/mol', 'Ni58':r'$\mu$mol/mol', 
            'Co59':r'$\mu$mol/mol', 'Cu63':r'$\mu$mol/mol', 'Zn64':r'$\mu$mol/mol', 
            'Rb85':r'$\mu$mol/mol', 'Sr88':'mmol/mol', 'Mo92':r'$\mu$mol/mol', 'Cd111':r'$\mu$mol/mol', 
            'Ba138':r'$\mu$mol/mol', 'Nd146':r'$\mu$mol/mol','U238':'nmol/mol'}
    return units[isotope]

#define paths
data_path=Path(os.getcwd())/"data"
fig_path=Path(os.getcwd())/"figures"


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


#plot the distribution of the isotopes through time
isotopes=['Li7', 'B11', 'Na23', 'Mg24', 'Al27', 'P31', 'S32', 'K39',
      'Mn55','Fe56','Ni58','Co59','Cu63', 'Zn64', 'Rb85', 'Sr88', 
      'Mo92','Cd111', 'Ba138', 'Nd146', 'U238']


#make Cenozoic epoch boundaries dict
epoch_boundaries={'Paleocene':(65.5, 55.8), 
                  'Eocene':(55.8, 33.9), 
                  'Oligocene':(33.9, 23.0), 
                  'Miocene':(23.0, 5.3), 
                  'Pliocene':(5.3, 1.8), 
                  'Pleistocene':(1.8,0.01)}

#order the foram species for colours Nuts first then Cibs. 
#units
"""""
for iso in isotopes:
    
    #set any negative data to np.nan
    negatives=foram_te_df[iso]<0
    foram_te_df.loc[negatives, iso]=np.nan
    
    #remove outliers
    out=find_outliers(foram_te_df[iso], mod=2.5)
    foram_te_df.loc[out, iso]=np.nan
    
    
    fig, ax =plt.subplots(figsize=(12, 6))
    p1=sns.scatterplot(data=foram_te_df, x='age_Ma', y=iso, 
                       hue='species_simple', style='core')
    #add epoch boundaries
    for epoch, bounds in epoch_boundaries.items():
        plt.axvline(bounds[0], alpha=0.2, color='grey', linestyle='--',)
    p1.legend(fontsize=6, ncol=2)
    
    plt.ylabel(f'{isotope_to_element(iso)}/Ca ({isotope_to_units(iso)})')
    plt.xlabel('Age (Ma)')
    fig.savefig(fig_path/f'{iso}_vs_age.png', dpi=300)
    print(f'{iso}: {sum(out)+sum(negatives)} outliers removed')
    

"""








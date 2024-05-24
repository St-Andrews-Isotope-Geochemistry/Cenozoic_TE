from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm


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

def make_stacked_plot(fig, ax, epoch_lines=True, adjust=0):

    for i, axis in enumerate(ax):
        #if even
        if i%2 == 0:
            axis.spines['right'].set_visible(False)
        else:
            axis.spines['left'].set_visible(False)
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')
        if i > 0:
            axis.spines['top'].set_visible(False)
        if i < len(ax)-1:
            axis.spines['bottom'].set_visible(False)
            plt.setp(axis.get_xticklabels(), visible=False)
            plt.setp(axis.get_xticklines(), visible=False)

        axis.set_facecolor('none')
        
        
        if epoch_lines:
            epoch_boundaries={'Paleocene':(65.5, 55.8), 
                  'Eocene':(55.8, 33.9), 
                  'Oligocene':(33.9, 23.0), 
                  'Miocene':(23.0, 5.3), 
                  'Pliocene':(5.3, 1.8), 
                  'Pleistocene':(1.8,0.01)}
            for epoch, bounds in epoch_boundaries.items():
                axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                            zorder=1)
        
    fig.subplots_adjust(hspace=adjust)


#make Cenozoic epoch boundaries dict
epoch_boundaries={'Paleocene':(65.5, 55.8), 
                  'Eocene':(55.8, 33.9), 
                  'Oligocene':(33.9, 23.0), 
                  'Miocene':(23.0, 5.3), 
                  'Pliocene':(5.3, 1.8), 
                  'Pleistocene':(1.8,0.01)}


#define paths
data_path=Path(os.getcwd())/"data"
fig_path=Path(os.getcwd())/"figures"
#import data
foram_df=pd.read_csv(data_path/"foram_dataframe.csv", index_col=0)
temp_df=pd.read_csv(data_path/"Meckler2022_temp.csv")
temp_df.dropna(subset=['temp_c'], inplace=True)
temp_df.sort_values('age_Ma', inplace=True)

Ca_df=pd.read_excel(data_path/"calcium_magnesium.xlsx", sheet_name='calcium')
Mg_df=pd.read_excel(data_path/"calcium_magnesium.xlsx", sheet_name='magnesium')
MgCasw_df=Ca_df.merge(Mg_df, on='age', suffixes=('_Ca', '_Mg'))
MgCasw_df['MgCa_sw']=MgCasw_df['median_Mg']/MgCasw_df['median_Ca']
cramer_temp=pd.read_csv(data_path/"Cramer_temp_d18O.csv")
cramer_MgCa=pd.read_csv(data_path/"Cramer_MgCa.csv")
Lear_MgCa=pd.read_csv(data_path/"Lear2015.csv")
#format data
foram_df['DeltaCO3']=foram_df['CO3']-foram_df['CO3']/foram_df['omega_logfit']
foram_df_Mg=foram_df.loc[pd.notnull(foram_df['Mg24'])]

foram_df_Mg['Mg/Ca_sw']=foram_df_Mg['Mg_sw']/foram_df_Mg['Ca_sw']

#reformat Lear
Lear_MgCa_form=Lear_MgCa[['Sample ID', 'Age (Ma)']].copy()
Lear_MgCa_form['Mg/Ca']=Lear_MgCa['O. umbonatus Mg/Ca'].copy()
Lear_MgCa_form['Mg/Ca_adj']=Lear_MgCa['O. umbonatus Mg/Ca_adj'].copy()
Lear_MgCa_form['Species']='Oridorsalis umbonatus'
Lear_MgCa_form.dropna(subset=['Mg/Ca'], inplace=True)

mund_df=Lear_MgCa[['Sample ID', 'Age (Ma)', 'C. mundulus Mg/Ca [mmol/mol]']].dropna(subset=['C. mundulus Mg/Ca [mmol/mol]'])
mund_df.rename(columns={'C. mundulus Mg/Ca [mmol/mol]':'Mg/Ca'}, inplace=True)
mund_df['Mg/Ca_adj']=mund_df['Mg/Ca'].copy()
mund_df['Species']='Cibicidoides mundulus'
Lear_MgCa_form=pd.concat([Lear_MgCa_form, mund_df], axis=0)

wuell_df=Lear_MgCa[['Sample ID', 'Age (Ma)', 'C. wuellerstorfi Mg/Ca [mmol/mol]']].dropna(subset=['C. wuellerstorfi Mg/Ca [mmol/mol]'])
wuell_df.rename(columns={'C. wuellerstorfi Mg/Ca [mmol/mol]':'Mg/Ca'}, inplace=True)
wuell_df['Mg/Ca_adj']=wuell_df['Mg/Ca'].copy()
wuell_df['Species']='Cibicidoides wuellerstorfi'
Lear_MgCa_form=pd.concat([Lear_MgCa_form, wuell_df], axis=0)




Lear_MgCa_form.loc[Lear_MgCa_form['Sample ID'].str.contains('690B'), 'ocean_basin']='Southern'
Lear_MgCa_form.loc[Lear_MgCa_form['Sample ID'].str.contains('806'), 'ocean_basin']='Pacific'
Lear_MgCa_form.loc[Lear_MgCa_form['Sample ID'].str.contains('926'), 'ocean_basin']='Atlantic'


cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('926'), 'ocean_basin']='Atlantic'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('806'), 'ocean_basin']='Pacific'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('690'), 'ocean_basin']='Southern'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('1088'), 'ocean_basin']='Atlantic'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('757'), 'ocean_basin']='Indian'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('747'), 'ocean_basin']='Southern'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('607'), 'ocean_basin']='Atlantic'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('573'), 'ocean_basin']='Pacific'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('761'), 'ocean_basin']='Indian'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('1171'), 'ocean_basin']='Southern'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('1218'), 'ocean_basin']='Pacific'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('689'), 'ocean_basin']='Southern'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('522'), 'ocean_basin']='Atlantic'
cramer_MgCa.loc[cramer_MgCa['Site'].str.contains('1209'), 'ocean_basin']='Pacific'


#species corrections
species_correct_Mg={'Cibicidoides havanensis':-0.4606, 'Cibicidoides grimsdalei':0.2872, 
                    'Cibicidoides laurisae':-0.5659, 'Cibicidoides subspiratus':0.4336, 
                    'Cibicidoides eocaenus': -0.1039, 
                    'Cibicidoides mundulus':0.0, 'Cibicidoides wuellerstorfi':-0.180458,
                    'Oridorsalis umbonatus':0.0}


foram_df_Mg=foram_df_Mg.loc[foram_df_Mg['species_simple'].isin(species_correct_Mg.keys())]

foram_df_Mg['Mg_species_corrected']=foram_df_Mg['Mg24']
for key, val in species_correct_Mg.items():
    foram_df_Mg.loc[foram_df_Mg['species_simple']==key, 'Mg_species_corrected']-=val

#what is the distribution of the species through time?
fig, ax =plt.subplots(figsize=(10, 6))
sns.stripplot(data=foram_df_Mg, x='age_Ma', hue='species_simple', 
              size=15, y='species_simple', ax=ax, alpha=.1, 
              jitter=False, linewidth=1, marker="s")
ax.set(ylabel=None)


cleaning_correct=0.909
foram_df_Mg['Mg_cleaning_corrected']=foram_df_Mg['Mg_species_corrected']*cleaning_correct


Mg_to_BWT={'Oridorsalis umbonatus':lambda Mg: (Mg-1.45)/0.1, 
            'Cibicidoides mundulus':lambda Mg: np.log(Mg/0.867)/0.109, 
            'Cibicidoides wuellerstorfi':lambda Mg: np.log(Mg/0.867)/0.109}





for species in pd.unique(foram_df_Mg['species_simple']):
    if species == 'Oridorsalis umbonatus':
        Mg_to_BWT = lambda Mg: (Mg-1.45)/0.1
    elif species == 'Cibicidoides mundulus': 
        Mg_to_BWT = lambda Mg: np.log(Mg/0.9)/0.11
    else:
        Mg_to_BWT = lambda Mg: np.log(Mg/0.867)/0.109
        
    foram_df_Mg.loc[foram_df_Mg['species_simple']==species, 'BWT_reduct']=Mg_to_BWT(foram_df_Mg.loc[foram_df_Mg['species_simple']==species, 'Mg_cleaning_corrected'])
    foram_df_Mg.loc[foram_df_Mg['species_simple']==species, 'BWT']=Mg_to_BWT(foram_df_Mg.loc[foram_df_Mg['species_simple']==species, 'Mg_species_corrected'])

    cramer_MgCa.loc[cramer_MgCa['Species']==species, 'BWT_reduct']=Mg_to_BWT(cramer_MgCa.loc[cramer_MgCa['Species']==species, 'Mg/Ca'])
    cramer_MgCa.loc[cramer_MgCa['Species']==species, 'BWT_reduct']=Mg_to_BWT(cramer_MgCa.loc[cramer_MgCa['Species']==species, 'Mg/Ca adj'])

    Lear_MgCa_form.loc[Lear_MgCa_form['Species']==species, 'BWT_reduct']=Mg_to_BWT(Lear_MgCa_form.loc[Lear_MgCa_form['Species']==species, 'Mg/Ca_adj'])



foram_df_BWT=foram_df_Mg.dropna(subset=['BWT'])

foram_df_BWT_Owullmund=foram_df_BWT.loc[foram_df_BWT['species_simple'].isin(['Cibicidoides mundulus', 'Cibicidoides wuellerstorfi', 'Oridorsalis umbonatus'])]
foram_df_BWT_other=foram_df_BWT.loc[~foram_df_BWT['species_simple'].isin(['Cibicidoides mundulus', 'Cibicidoides wuellerstorfi', 'Oridorsalis umbonatus'])]

## plot with Meckler data showing the correction
fig, ax =plt.subplots(figsize=(8, 16), nrows=4, sharex=True)
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='s',
                label='Meckler et al. (2022) clumped isotopes', ax=ax[0], legend=False)
ax[0].set_ylabel('Clumped isotopes ($^{\circ}$C)')
counter=1
ax[0].set_ylim(-3, 21)
ax[0].set_yticks([-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
for key in pd.unique(foram_df_BWT_Owullmund['species_simple']):
    p2=sns.lineplot(data=foram_df_BWT_Owullmund.loc[foram_df_BWT_Owullmund['species_simple']==key], x='age_Ma', y='BWT', 
                    marker='o', ax=ax[counter], hue='ocean_basin')
    #p2=sns.lineplot(data=foram_df_BWT.loc[foram_df_BWT['species_simple']==key], x='age_Ma', y='BWT_reduct', 
    #                marker='o', ax=ax[counter], linestyle='--', alpha=.5, 
    #                hue='ocean_basin')
    ax[counter].set_ylabel(key+' ($^{\circ}$C)')
    ax[counter].legend(loc='lower left')
    ax[counter].set_ylim(-3, 21)
    if key == 'Cibicidoides wuellerstorfi':
        ax[counter].set_yticks([-2, 0, 2, 4, 6])
    elif key == 'Cibicidoides mundulus':
        ax[counter].set_yticks([0, 2, 4, 6, 8, 10])
    else:
        ax[counter].set_yticks([-2, 0, 2, 4, 6, 8, 10, 12, 14, 16])

    counter+=1
    

ax[-1].set_xlim(0, 65)
ax[-1].invert_xaxis()
ax[-1].set_xlabel('Age (Ma)')
plt.subplots_adjust(hspace=-0.7)
make_stacked_plot(ax, epoch_lines=True)
fig.savefig(fig_path/'BWT_1.png', dpi=300)



## plot with Meckler data on same axis
hue_order=['Atlantic', 'Pacific', 'Southern']
fig, ax =plt.subplots(figsize=(10, 8))
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^',
                label='Meckler et al. (2022) clumped isotopes', mfc='white', mec='black')
p4=sns.lineplot(data=cramer_temp, x='Age', y='Temperature', color='firebrick', 
                label=r'Cramer et al. (2016) $\delta^{18}$O', ax=ax)
p2=sns.lineplot(data=foram_df_BWT_Owullmund, 
                x='age_Ma', y='BWT', hue='ocean_basin', palette='Set2', markers=True, 
                hue_order=hue_order, style='species_simple')
p3=sns.scatterplot(data=foram_df_BWT_other, x='age_Ma', y='BWT', hue='ocean_basin', palette='Set2', 
                   hue_order=hue_order, marker='P', legend=False)
# manually add legend entry for other species using only the marker (not the hue)
ax.plot([], [], marker='P', color='black', label='Other (corrected to C. mundulus)', ls='None')


ax.set_xlim(0, 65)
ax.legend()
ax.invert_xaxis()

plt.title('Benthic Mg/Ca')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlabel('Age (Ma)')                
fig.savefig(fig_path/'BWT_2.png', dpi=300)

## plot with Meckler data on same axis using reductive correction
hue_order=['Atlantic', 'Pacific', 'Southern']
fig, ax =plt.subplots(figsize=(10, 8))
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^',
                label='Meckler et al. (2022) clumped isotopes', mfc='white', mec='black')
p4=sns.lineplot(data=cramer_temp, x='Age', y='Temperature', color='grey', 
                label=r'Cramer et al. (2016) $\delta^{18}$O', ax=ax)
p2=sns.lineplot(data=foram_df_BWT_Owullmund, 
                x='age_Ma', y='BWT_reduct', hue='ocean_basin', palette='Set2', markers=True, 
                hue_order=hue_order, style='species_simple')
p3=sns.scatterplot(data=foram_df_BWT_other, x='age_Ma', y='BWT_reduct', hue='ocean_basin', palette='Set2', 
                   hue_order=hue_order, marker='P', legend=False)
# manually add legend entry for other species using only the marker (not the hue)
ax.plot([], [], marker='P', color='black', label='Other (corrected to C. mundulus)', ls='None')

ax.legend()
ax.set_xlim(0, 65)
ax.invert_xaxis()

plt.title('Benthic Mg/Ca with reductive cleaning correction')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlabel('Age (Ma)')                
fig.savefig(fig_path/'BWT_2_reductive.png', dpi=300)











## plot with Meckler data on same axis with Mg_sw
hue_order=['Atlantic', 'Pacific', 'Southern']
fig, ax =plt.subplots(figsize=(10, 10), nrows=2, sharex=True, height_ratios=[2, 1])
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^',
                label='Meckler et al. (2022) clumped isotopes', mfc='white', mec='black', ax=ax[0])
p4=sns.lineplot(data=cramer_temp, x='Age', y='Temperature', color='firebrick', 
                label=r'Cramer et al. (2016) $\delta^{18}$O', ax=ax[0])
p2=sns.lineplot(data=foram_df_BWT_Owullmund, 
                x='age_Ma', y='BWT', hue='ocean_basin', palette='Set2', markers=True, 
                hue_order=hue_order, style='species_simple', ax=ax[0])
p3=sns.scatterplot(data=foram_df_BWT_other, x='age_Ma', y='BWT', hue='ocean_basin', palette='Set2', 
                   hue_order=hue_order, marker='P', legend=False, ax=ax[0])
# manually add legend entry for other species using only the marker (not the hue)
ax[0].plot([], [], marker='P', color='black', label='Other (corrected to C. mundulus)', ls='None')
ax[0].legend(loc='upper right', fontsize=8, ncol=2)
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[0].set_ylim(-3, 21)   

sns.lineplot(data=foram_df_Mg, x='age_Ma', y='Mg/Ca_sw', ax=ax[1], color='black', label='Mg/Ca_sw', 
             legend=False)
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/kg)') 

ax[1].set_xlabel('Age (Ma)')     
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
plt.subplots_adjust(hspace=-0.3)
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)
fig.savefig(fig_path/'BWT_Mg_sw.png', dpi=300)






fig, ax =plt.subplots(figsize=(10, 8))
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^',
                label='Meckler et al. (2022) clumped isotopes', 
                mfc='white', mec='black', ax=ax)
sns.scatterplot(data=Lear_MgCa_form, x='Age (Ma)', y='BWT_reduct', 
                hue='ocean_basin', style='Species', palette='Set2', ax=ax)
ax.set_xlim(0, 65)
ax.invert_xaxis()
ax.set_xlabel('Age (Ma)')
ax.set_ylabel('BWT ($^{\circ}$C)')
plt.title('Lear data')
fig.savefig(fig_path/'Lear_BWT.png', dpi=300)


Lear_MgCa_form_Atl=Lear_MgCa_form.loc[Lear_MgCa_form['ocean_basin']=='Atlantic']

fig, ax =plt.subplots(figsize=(10, 8))
p1=sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^',
                label='Meckler et al. (2022) clumped isotopes',
                mfc='white', mec='black', ax=ax)
sns.scatterplot(data=cramer_MgCa, x='Age', y='BWT_reduct', 
                hue='ocean_basin', style='Species', palette='Set2', ax=ax)
ax.set_xlim(0, 65)
ax.invert_xaxis()
ax.set_xlabel('Age (Ma)')
ax.set_ylabel('BWT ($^{\circ}$C)')
plt.title('Cramer data')
fig.savefig(fig_path/'Cramer_BWT.png', dpi=300)



cramer_MgCa_Atl=cramer_MgCa.loc[cramer_MgCa['ocean_basin']=='Atlantic']


## smooth temp data
# Assuming temp_df is your DataFrame and it has columns 'age_Ma', 'temp_c', and 'err'
# temp_df = pd.DataFrame({'age_Ma': ..., 'temp_c': ..., 'err': ...})

# Create a loess fit
lowess = sm.nonparametric.lowess(temp_df['temp_c'], temp_df['age_Ma'], frac=0.2, it=3, delta=0.0, is_sorted=False, missing='drop', return_sorted=True)
# Create a new DataFrame with the results
lowess_df = pd.DataFrame(lowess, columns=['age_Ma', 'temp_c_smooth'])
# Merge the original data with the lowess smoothed data
merged_df = pd.merge(temp_df, lowess_df, on='age_Ma')
# Plot the original data
plt.errorbar(merged_df['age_Ma'], merged_df['temp_c'], yerr=merged_df['err'], linestyle='None', marker='o')
# Plot the smoothed data
plt.plot(merged_df['age_Ma'], merged_df['temp_c_smooth'], color='red')
plt.show()


#interpolate onto foram data
foram_df_Mg['temp_c_smooth']=np.nan
foram_df_Mg['temp_c_smooth']=np.interp(foram_df_Mg['age_Ma'], temp_df['age_Ma'], lowess_df['temp_c_smooth'])


foram_df_Mg_cibs=foram_df_Mg.loc[foram_df_Mg['species_simple'].str.contains('Cibicidoides')]
fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_Mg_cibs, x='age_Ma', y='Mg_species_corrected', hue='ocean_basin', ax=ax)



## fit Lear models to Cibs linear

from scipy.optimize import curve_fit

def expo_fit(X, A, H, B):
    MgCa_sw, BWT = X
    return A*MgCa_sw**H*np.exp(B*BWT)

def expo_fit_predict(X, A, H, B):
    MgCa, MgCa_sw = X
    return np.log(MgCa/(A*MgCa_sw**H))/B


def lin_fit(X, c, m, H):
    MgCa_sw, BWT = X
    return (c+m*BWT)*MgCa_sw**H

def lin_fit_predict(X, c, m, H):
    MgCa, MgCa_sw = X
    return (MgCa/MgCa_sw**H-c)/m




xdata=(foram_df_Mg_cibs['Mg/Ca_sw'], foram_df_Mg_cibs['temp_c_smooth'].values)
ydata=foram_df_Mg_cibs['Mg_species_corrected'].values

parameters, covariance = curve_fit(lin_fit, xdata, ydata)

xdata_predict=(foram_df_Mg_cibs['Mg_species_corrected'], foram_df_Mg_cibs['Mg/Ca_sw'])
foram_df_Mg_cibs['BWT']=lin_fit_predict(xdata_predict, *parameters)

x_predict=np.linspace(foram_df_Mg_cibs['Mg_species_corrected'].min(), foram_df_Mg_cibs['Mg_species_corrected'].max())
ypredict_2=lin_fit_predict((x_predict, np.array([2]*len(x_predict))), *parameters)
ypredict_3=lin_fit_predict((x_predict, np.array([3]*len(x_predict))), *parameters)
ypredict_4=lin_fit_predict((x_predict, np.array([4]*len(x_predict))), *parameters)
ypredict_5=lin_fit_predict((x_predict, np.array([5]*len(x_predict))), *parameters)


fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_Mg_cibs, y='temp_c_smooth', x='Mg_species_corrected', hue='Mg/Ca_sw', ax=ax)
plt.plot(x_predict, ypredict_2, label='Mg/Ca_sw=2')
plt.plot(x_predict, ypredict_3 ,label='Mg/Ca_sw=3')
plt.plot(x_predict, ypredict_4, label='Mg/Ca_sw=4')
plt.plot(x_predict, ypredict_5,label='Mg/Ca_sw=5')
plt.legend()

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Mg_cibs, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Cibs only. Lear et al., 2015, linear')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)


parameters_dict={'Lear_linear_Cibs':parameters}





## fit Lear models to Cibs exponential

parameters, covariance = curve_fit(expo_fit, xdata, ydata)

foram_df_Mg_cibs['BWT']=expo_fit_predict(xdata_predict, *parameters)

ypredict_2=expo_fit_predict((x_predict, np.array([2]*len(x_predict))), *parameters)
ypredict_3=expo_fit_predict((x_predict, np.array([3]*len(x_predict))), *parameters)
ypredict_4=expo_fit_predict((x_predict, np.array([4]*len(x_predict))), *parameters)
ypredict_5=expo_fit_predict((x_predict, np.array([5]*len(x_predict))), *parameters)


fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_Mg_cibs, y='temp_c_smooth', x='Mg_species_corrected', hue='Mg/Ca_sw', ax=ax)
plt.plot(x_predict, ypredict_2, label='Mg/Ca_sw=2')
plt.plot(x_predict, ypredict_3 ,label='Mg/Ca_sw=3')
plt.plot(x_predict, ypredict_4, label='Mg/Ca_sw=4')
plt.plot(x_predict, ypredict_5,label='Mg/Ca_sw=5')
plt.legend()

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Mg_cibs, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Cibs only. Lear et al., 2015, exponential')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)


parameters_dict['Lear_exponential_Cibs']= parameters











##Add Lear data
Lear_MgCa_form_Cibs=Lear_MgCa_form.loc[Lear_MgCa_form['Species'].str.contains('Cibicidoides')]
Lear_concat_df=pd.DataFrame()


Lear_concat_df['core']=Lear_MgCa_form['Sample ID']
Lear_concat_df['core']=np.where(Lear_concat_df['core'].str.contains('690B'), '690B', Lear_concat_df['core'])
Lear_concat_df['core']=np.where(Lear_concat_df['core'].str.contains('806B'), '806B', Lear_concat_df['core'])
Lear_concat_df['age_Ma']=Lear_MgCa_form['Age (Ma)']
Lear_concat_df['Mg24']=Lear_MgCa_form['Mg/Ca']
Lear_concat_df['Mg_species_corrected']=Lear_MgCa_form['Mg/Ca']
Lear_concat_df['species_simple']=Lear_MgCa_form['Species']
Lear_concat_df['ocean_basin']=Lear_MgCa_form['ocean_basin']
Lear_concat_df['source']='Lear et al. (2015)'
Lear_concat_df['temp_c_smooth']=np.nan  
Lear_concat_df['temp_c_smooth']=np.interp(Lear_concat_df['age_Ma'], temp_df['age_Ma'], lowess_df['temp_c_smooth'])
Lear_concat_df['Mg_sw']=np.interp(Lear_concat_df['age_Ma'], Mg_df['age'], Mg_df['median'])
Lear_concat_df['Ca_sw']=np.interp(Lear_concat_df['age_Ma'], Ca_df['age'], Ca_df['median'])
Lear_concat_df['Mg/Ca_sw']=Lear_concat_df['Mg_sw']/Lear_concat_df['Ca_sw']

foram_df_wLear=pd.concat([foram_df_Mg, Lear_concat_df], axis=0)
foram_df_wLear_cibs=foram_df_wLear.loc[foram_df_wLear['species_simple'].str.contains('Cibicidoides')]
foram_df_wLear_cibs_Atl=foram_df_wLear_cibs.loc[foram_df_wLear_cibs['ocean_basin']=='Atlantic']

hue_order=['Atlantic', 'Pacific', 'Southern']
fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_wLear_cibs, x='age_Ma', y='Mg_species_corrected', ax=ax,
                palette='pastel', hue='ocean_basin', hue_order=hue_order)
sns.scatterplot(data=Lear_concat_df, x='age_Ma', y='Mg_species_corrected', ax=ax, 
                hue='ocean_basin', hue_order=hue_order)
ax.set_xlim(0, 65)
ax.invert_xaxis()




## fit with +Lear data (Atl Cibs only)

#fit to only Atlantic
xdata=(foram_df_wLear_cibs_Atl['Mg/Ca_sw'], foram_df_wLear_cibs_Atl['temp_c_smooth'].values)
ydata=foram_df_wLear_cibs_Atl['Mg_species_corrected'].values

parameters, covariance = curve_fit(lin_fit, xdata, ydata)
#predict with all basins
xdata_predict=(foram_df_wLear_cibs['Mg_species_corrected'], foram_df_wLear_cibs['Mg/Ca_sw'])
foram_df_wLear_cibs['BWT']=lin_fit_predict(xdata_predict, *parameters)

x_predict=np.linspace(foram_df_wLear_cibs_Atl['Mg_species_corrected'].min(), foram_df_wLear_cibs_Atl['Mg_species_corrected'].max())
ypredict_2=lin_fit_predict((x_predict, np.array([2]*len(x_predict))), *parameters)
ypredict_3=lin_fit_predict((x_predict, np.array([3]*len(x_predict))), *parameters)
ypredict_4=lin_fit_predict((x_predict, np.array([4]*len(x_predict))), *parameters)
ypredict_5=lin_fit_predict((x_predict, np.array([5]*len(x_predict))), *parameters)


fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_wLear_cibs_Atl, y='temp_c_smooth', x='Mg_species_corrected', hue='Mg/Ca_sw', ax=ax)
plt.plot(x_predict, ypredict_2, label='Mg/Ca_sw=2')
plt.plot(x_predict, ypredict_3 ,label='Mg/Ca_sw=3')
plt.plot(x_predict, ypredict_4, label='Mg/Ca_sw=4')
plt.plot(x_predict, ypredict_5,label='Mg/Ca_sw=5')
plt.legend()

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_wLear_cibs, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Cibs only. Lear et al., 2015, linear, w Lear data')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)


parameters_dict['Lear_linear_Cibs_wLearData']=parameters



## Use simple lin regress instead on Cibs

xdata=foram_df_wLear_cibs_Atl[['Mg/Ca_sw', 'Mg_species_corrected']].values
ydata=foram_df_wLear_cibs_Atl['temp_c_smooth'].values
xdata=sm.add_constant(xdata)
model = sm.OLS(ydata, xdata)
results=model.fit()
xdata_predict=foram_df_wLear_cibs[['Mg/Ca_sw', 'Mg_species_corrected']].values
xdata_predict=sm.add_constant(xdata_predict)
predictions=results.predict(xdata_predict)
parameters=results.params
foram_df_wLear_cibs_Atl['BWT']=predictions

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_wLear_cibs, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Cibs only. OLS, w Lear data')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)

parameters_dict['OLS_Cibs_wLearData']=parameters


## OLS with Li7 (no Lear data)


fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Mg_cibs, x='age_Ma', y='BWT', hue='Li7', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
#fig.suptitle('Cibs only. OLS, w Lear data')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)




xdata=foram_df_Mg_cibs[['Mg/Ca_sw', 'Mg_species_corrected', 'Li7']].values
ydata=foram_df_Mg_cibs['temp_c_smooth'].values
xdata=sm.add_constant(xdata)
model = sm.OLS(ydata, xdata)
results=model.fit()
xdata_predict=foram_df_Mg_cibs[['Mg/Ca_sw', 'Mg_species_corrected', 'Li7']].values
xdata_predict=sm.add_constant(xdata)
predictions=results.predict(xdata_predict)
parameters=results.params
foram_df_Mg_cibs['BWT']=predictions

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Mg_cibs, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Cibs only. OLS, w Li7')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)






## Oridorsalis umbonatus calibration


foram_df_Orid=foram_df_Mg.loc[foram_df_Mg['species_simple']=='Oridorsalis umbonatus']
fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_Orid, x='age_Ma', y='Mg_species_corrected', ax=ax, hue='ocean_basin')
ax.set_xlim(0, 65)
ax.invert_xaxis()







xdata=(foram_df_Orid['Mg/Ca_sw'], foram_df_Orid['temp_c_smooth'].values)
ydata=foram_df_Orid['Mg_species_corrected'].values

parameters, covariance = curve_fit(lin_fit, xdata, ydata)

xdata_predict=(foram_df_Orid['Mg_species_corrected'], foram_df_Orid['Mg/Ca_sw'])
foram_df_Orid['BWT']=lin_fit_predict(xdata_predict, *parameters)

x_predict=np.linspace(foram_df_Orid['Mg_species_corrected'].min(), foram_df_Orid['Mg_species_corrected'].max())
ypredict_2=lin_fit_predict((x_predict, np.array([2]*len(x_predict))), *parameters)
ypredict_3=lin_fit_predict((x_predict, np.array([3]*len(x_predict))), *parameters)
ypredict_4=lin_fit_predict((x_predict, np.array([4]*len(x_predict))), *parameters)
ypredict_5=lin_fit_predict((x_predict, np.array([5]*len(x_predict))), *parameters)


fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_Orid, y='temp_c_smooth', x='Mg_species_corrected', hue='Mg/Ca_sw', ax=ax)
plt.plot(x_predict, ypredict_2, label='Mg/Ca_sw=2')
plt.plot(x_predict, ypredict_3 ,label='Mg/Ca_sw=3')
plt.plot(x_predict, ypredict_4, label='Mg/Ca_sw=4')
plt.plot(x_predict, ypredict_5,label='Mg/Ca_sw=5')
plt.legend()

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Orid, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Oridisalis only. Lear et al., 2015, linear')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)


parameters_dict={'Lear_linear_Oridisalis':parameters}

## With Lear data

foram_df_wLear_Orid=foram_df_wLear.loc[foram_df_wLear['species_simple']=='Oridorsalis umbonatus']

xdata=(foram_df_wLear_Orid['Mg/Ca_sw'], foram_df_wLear_Orid['temp_c_smooth'].values)
ydata=foram_df_wLear_Orid['Mg_species_corrected'].values

parameters, covariance = curve_fit(lin_fit, xdata, ydata)

xdata_predict=(foram_df_wLear_Orid['Mg_species_corrected'], foram_df_wLear_Orid['Mg/Ca_sw'])
foram_df_wLear_Orid['BWT']=lin_fit_predict(xdata_predict, *parameters)

x_predict=np.linspace(foram_df_wLear_Orid['Mg_species_corrected'].min(), foram_df_wLear_Orid['Mg_species_corrected'].max())
ypredict_2=lin_fit_predict((x_predict, np.array([2]*len(x_predict))), *parameters)
ypredict_3=lin_fit_predict((x_predict, np.array([3]*len(x_predict))), *parameters)
ypredict_4=lin_fit_predict((x_predict, np.array([4]*len(x_predict))), *parameters)
ypredict_5=lin_fit_predict((x_predict, np.array([5]*len(x_predict))), *parameters)


fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_wLear_Orid, y='temp_c_smooth', x='Mg_species_corrected', hue='Mg/Ca_sw', ax=ax)
plt.plot(x_predict, ypredict_2, label='Mg/Ca_sw=2')
plt.plot(x_predict, ypredict_3 ,label='Mg/Ca_sw=3')
plt.plot(x_predict, ypredict_4, label='Mg/Ca_sw=4')
plt.plot(x_predict, ypredict_5,label='Mg/Ca_sw=5')
plt.legend()

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_wLear_Orid, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Oridisalis only. Lear et al., 2015, linear, w Lear data')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)


parameters_dict={'Lear_linear_Oridisalis_wLearData':parameters}


## OLS on Oridisalis



xdata=foram_df_wLear_Orid[['Mg/Ca_sw', 'Mg_species_corrected']].values
ydata=foram_df_wLear_Orid['temp_c_smooth'].values
xdata=sm.add_constant(xdata)
model = sm.OLS(ydata, xdata)
results=model.fit()
xdata_predict=foram_df_wLear_Orid[['Mg/Ca_sw', 'Mg_species_corrected']].values
xdata_predict=sm.add_constant(xdata_predict)
predictions=results.predict(xdata_predict)
parameters=results.params
foram_df_wLear_Orid['BWT']=predictions

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_wLear_Orid, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Oridisalis only. OLS, w Lear data')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)

parameters_dict['OLS_Oridisalis_wLearData']=parameters


## Oridisalis w Li7



xdata=foram_df_Orid[['Mg/Ca_sw', 'Mg_species_corrected', 'Li7']].values
ydata=foram_df_Orid['temp_c_smooth'].values
xdata=sm.add_constant(xdata)
model = sm.OLS(ydata, xdata)
results=model.fit()
xdata_predict=foram_df_Orid[['Mg/Ca_sw', 'Mg_species_corrected', 'Li7']].values
xdata_predict=sm.add_constant(xdata_predict)
predictions=results.predict(xdata_predict)
parameters=results.params
foram_df_Orid['BWT']=predictions

fig, ax =plt.subplots(figsize=(10, 8), nrows=2, sharex=True, height_ratios=[2, 1])
sns.scatterplot(data=foram_df_Orid, x='age_Ma', y='BWT', hue='ocean_basin', ax=ax[0])
sns.lineplot(data=temp_df, x='age_Ma', y='temp_c', color='black', marker='^', ax=ax[0])
sns.lineplot(data=merged_df, x='age_Ma', y='temp_c_smooth', color='black', 
             ax=ax[0], linestyle='--')
sns.lineplot(data=MgCasw_df, x='age', y='MgCa_sw', color='red', ax=ax[1])
ax[0].set_ylabel('BWT ($^{\circ}$C)')
ax[1].set_ylabel(r'Mg/Ca$_{sw}$ ($\mu$mol/mol)')
ax[1].set_xlabel('Age (Ma)')
ax[1].set_xlim(0, 65)
ax[1].invert_xaxis()
fig.suptitle('Oridisalis only. OLS, w Li7')
make_stacked_plot(fig, ax, epoch_lines=True, adjust=-0.3)

parameters_dict['OLS_Oridisalis_wLi7']=parameters



## add cramer data

#remove 926, 690 and 806 to prevent overlap with Lear
Lear_cores=['ODP0926', 'ODP0690', 'ODP0806']
cramer_MgCa_noLear=cramer_MgCa.loc[~cramer_MgCa['Site'].isin(Lear_cores)]  

cramer_MgCa_noLear.rename(columns={'Age':'age_Ma', 'Mg/Ca':'Mg24', 
                                   'Species':'species_simple', 'Site':'core'}, inplace=True)

cramer_MgCa_noLear['Mg_species_corrected']=cramer_MgCa_noLear['Mg24']

cramer_MgCa_noLear_umbwuelmun=cramer_MgCa_noLear.loc[cramer_MgCa_noLear['species_simple'].isin(['Cibicidoides mundulus', 'Cibicidoides wuellerstorfi', 'Oridorsalis umbonatus'])]

cramer_MgCa_noLear_umbwuelmun['source']='Cramer et al. (2011)'

foram_df_wLear_wcramer=pd.concat([foram_df_wLear, cramer_MgCa_noLear_umbwuelmun], axis=0)




fig, ax =plt.subplots(figsize=(10, 8))
sns.scatterplot(data=foram_df_wLear_wcramer, x='age_Ma', y='Mg_species_corrected', 
                hue='species_simple', ax=ax, style='source')
plt.legend(ncol=2)





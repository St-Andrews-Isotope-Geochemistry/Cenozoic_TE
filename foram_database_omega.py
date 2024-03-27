
from pathlib import Path
import os
import pandas as pd
import kgen
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm

def depth_m_to_pressure_bar(depth):
    pressure_bar=1023.6*9.80665*depth*10**-5
    return pressure_bar
def linregress(df, x_transform=None , y_transform=None, features='B_Ca_umolmol', target='omega_c'):
    
    X=df[features].values
    y=df[target].values
    if x_transform != None:
        X=x_transform(X)
    if y_transform != None:
        y=y_transform(y)
    
    #Reshape and add intercept to bracketed values
    Xwint=np.empty(shape=(len(X), 2), dtype=np.float64)
    Xwint[:,0]=1
    Xwint[:, 1]=X
    #reshape true values
    Y=y.reshape(-1, 1)

    mdl = sm.OLS(Y, Xwint)
    res_ols = mdl.fit()

    return res_ols
def linpred(xpredict, mdl, x_transform=None, y_transform=None):
    
    
    if x_transform != None:
        xpredict=x_transform(xpredict)
    
    if type(xpredict)==float:
        xpredwint=np.array([1, xpredict])
    else:
        xpredwint=np.empty(shape=(len(xpredict), 2), dtype=np.float64)
        xpredwint[:,0]=1
        xpredwint[:, 1]=xpredict
    ypred=mdl.predict(xpredwint)
    
    if y_transform != None:
        ypred=y_transform(ypred)
    
    return ypred
def omega_depth_correction(omega, depth, temp=2, newdepth=0):
    pressure_coefficients = kgen.coefs.K_presscorr_coefs
    pressure_correction = kgen.K_functions.calc_pressure_correction(coefficients=pressure_coefficients['KspC'], p_bar=depth_m_to_pressure_bar(depth), temp_c=temp)
    omega0=omega*pressure_correction
    pressure_correction = kgen.K_functions.calc_pressure_correction(coefficients=pressure_coefficients['KspC'], p_bar=depth_m_to_pressure_bar(newdepth), temp_c=temp)
    return omega0/pressure_correction


#Read Yu and Elderfield data
data_path=Path(os.getcwd())/"data"
file_path=data_path/"Yu_Elderfield_calibration.xlsx"
Yu_Elder_df=pd.read_excel(file_path)

#convert depth to pressure
Yu_Elder_df['pressure_bar']=depth_m_to_pressure_bar(Yu_Elder_df['water_depth_m'].values)



#set Ca and Mg concentrations and salinity
calcium_conc=0.01
mag_conc=0.053
salinity=35

#calculate KspC and KspA from kgen
dict_ksp=kgen.calc_Ks(["KspC","KspA"], temp_c=Yu_Elder_df['temp_c'].values, sal=salinity, p_bar=Yu_Elder_df['pressure_bar'].values, 
                       calcium=calcium_conc, magnesium=mag_conc)

Yu_Elder_df['KspC']=dict_ksp['KspC']
Yu_Elder_df['KspA']=dict_ksp['KspA']

#calculate CO3sat_c and CO3sat_c
Yu_Elder_df['CO3sat_c']=Yu_Elder_df['KspC']/calcium_conc*10**6
Yu_Elder_df['CO3sat_a']=Yu_Elder_df['KspA']/calcium_conc*10**6

#calculate omega
Yu_Elder_df['omega_c']=(Yu_Elder_df['DCO3_c']+Yu_Elder_df['CO3sat_c'])/Yu_Elder_df['CO3sat_c']
Yu_Elder_df['omega_a']=(Yu_Elder_df['DCO3_a']+Yu_Elder_df['CO3sat_a'])/Yu_Elder_df['CO3sat_a']



#plot B/Ca vs omega_c
fig, ax = plt.subplots()
sns.scatterplot(data=Yu_Elder_df, x='B_Ca_umolmol', y='omega_c', hue='species')
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
plt.show()
Yu_Elder_df['log(omega_c)']=np.log(Yu_Elder_df['omega_c'])


#plot side by side with log transformed omega_c
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
sns.scatterplot(ax=ax[0], data=Yu_Elder_df, x='B_Ca_umolmol', y='omega_c', hue='species')
ax[0].set_title('raw')
ax[0].set_xlabel(r'B/Ca ($\mu$mol/mol)')
sns.scatterplot(ax=ax[1], data=Yu_Elder_df, x='B_Ca_umolmol', y='log(omega_c)', hue='species')
ax[1].set_title('transformed')
ax[1].set_xlabel(r'B/Ca ($\mu$mol/mol)')


#get species names
species_names=pd.unique(Yu_Elder_df['species'])

#make dictionaries to store linear and log fits
species_fit_dict={}

#Make colormap for species
colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][0:len(species_names)]
species_col=dict(zip(species_names, colors))

#cycle through species and fit linear and log fits, and plot
fig, ax =plt.subplots(figsize=(6, 6))
for species in species_names:
    df=Yu_Elder_df.loc[Yu_Elder_df['species']==species]
    
    X=df['B_Ca_umolmol']
    Y=df['omega_c']
    
    xpredict=np.linspace(min(X), max(X), 30)
    
    linmdl= linregress(df)
    logmdl = linregress(df, y_transform=np.log)
    
    species_fit_dict[species]={'lin':linmdl, 'log':logmdl}
    
    species_short=species.split(' ')[1]
    
    ax.scatter(X, Y,color=species_col[species], label=species)
    ax.plot(xpredict, linpred(xpredict, linmdl), color=species_col[species], 
               label=f'r_sq = {linmdl.rsquared:.2f}')
    ax.plot(xpredict, linpred(xpredict, logmdl, y_transform=np.exp), color=species_col[species], 
               label=f'r_sq = {logmdl.rsquared:.2f}', ls='dashed')
    
ax.legend()
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
ax.set_ylabel(r'$\Omega_c$')






## Foram data   
#import foram data
excel_file_path = data_path/"foram_dataset.xlsx"

# Read all sheets from the Excel file into a dictionary of DataFrames
excel_data = pd.read_excel(excel_file_path, sheet_name=None)

foram_df=pd.DataFrame()
# Iterate through the sheets and store each DataFrame in the dictionary
for sheet_name, core_df in excel_data.items():
    core_df.insert(0, 'core', sheet_name)
    foram_df=pd.concat([foram_df, core_df])

foram_df.reset_index(drop=True, inplace=True)


#1262 paleodepth
foram_df.loc[foram_df['core']=='1262', 'palaeo_depth_m']=np.interp(foram_df.loc[foram_df['core']=='1262', 'age_Ma'], [56, 66], [3000, 3500])

#drop 690B
foram_df.drop(index=foram_df.index[foram_df['core']=='690B'], inplace=True)

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




foram_df['omega_linfit']=np.nan
foram_df['omega_logfit']=np.nan

for species in species_names:
    
    #perform linear and log predictions
    foram_df['omega_linfit']=np.where(foram_df['species_simple']==species, 
                                      linpred(foram_df['B11'].values, species_fit_dict[species]['lin']), 
                                      foram_df['omega_linfit'])
    
    foram_df['omega_logfit']=np.where(foram_df['species_simple']==species,
                                        linpred(foram_df['B11'].values, species_fit_dict[species]['log'], y_transform=np.exp),
                                        foram_df['omega_logfit'])
    


#isolate Cibicidoides for plotting
cibs_df=foram_df.dropna(subset='species', axis=0)
cibs_df=cibs_df.loc[cibs_df['species_simple'].isin(species_names)]

#plot
fig, ax =plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex=True)
sns.scatterplot(ax=ax[0], data=cibs_df, x='age_Ma', y='omega_logfit', hue='core', style='species_simple', legend=False)
sns.scatterplot(ax=ax[1], data=cibs_df, x='age_Ma', y='B11', hue='core', style='species_simple')
plt.legend(fontsize=6)






#do depth correction
foram_df['omega_linfit_0']=omega_depth_correction(foram_df['omega_linfit'], foram_df['palaeo_depth_m'])
foram_df['omega_logfit_0']=omega_depth_correction(foram_df['omega_logfit'], foram_df['palaeo_depth_m'])
foram_df['omega_linfit_3000']=omega_depth_correction(foram_df['omega_linfit'], foram_df['palaeo_depth_m'], newdepth=3000)
foram_df['omega_logfit_3000']=omega_depth_correction(foram_df['omega_logfit'], foram_df['palaeo_depth_m'], newdepth=3000)


#save to disk
foram_df.to_csv(data_path/"foram_dataframe.csv")


#isolate Cibicidoides for plotting
cibs_df=foram_df.dropna(subset='species', axis=0)
cibs_df=cibs_df.loc[cibs_df['species_simple'].isin(species_names)]

#plot
fig, ax =plt.subplots(nrows=2, ncols=1,  figsize=(12, 8), sharex=True, sharey=True)
sns.scatterplot(ax=ax[0], data=cibs_df, x='age_Ma', y='omega_logfit', hue='core')
sns.scatterplot(ax=ax[1], data=cibs_df, x='age_Ma', y='omega_logfit_3000', hue='core')
plt.show()






## Ca and Mg

#read in Ca and Mg data
Ca_df=pd.read_excel(data_path/"calcium_magnesium.xlsx", sheet_name='calcium')
Mg_df=pd.read_excel(data_path/"calcium_magnesium.xlsx", sheet_name='magnesium')


#perform linear interpolation on age_Ma in foram_df to get calcium and magnesium values
foram_df['Ca_sw']=np.interp(foram_df['age_Ma'], Ca_df['age'], Ca_df['median'])
foram_df['Mg_sw']=np.interp(foram_df['age_Ma'], Mg_df['age'], Mg_df['median'])

#plot
fig, ax =plt.subplots(figsize=(10, 6))
sns.scatterplot(data=foram_df, x='age_Ma', y='Ca_sw', color='black')
sns.scatterplot(data=foram_df, x='age_Ma', y='Mg_sw', color='red')
sns.lineplot(data=Ca_df, x='age', y='median', color='black')
sns.lineplot(data=Mg_df, x='age', y='median', color='red')
ax.set_ylabel(r'Concentration ($\mu$mol/kg)')
ax.set_xlabel('Age (Ma)')
ax.set_xlim(0, 70)
plt.legend(['Calcium', 'Magnesium'])
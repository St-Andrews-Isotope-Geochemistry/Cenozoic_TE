
from pathlib import Path
import os
import pandas as pd
import kgen
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import time

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
fig_path=Path(os.getcwd())/"figures"
file_path=data_path/"Yu_Elderfield_calibration.xlsx"
Yu_Elder_df=pd.read_excel(file_path)
Yu_Elder_df['reference']='Yu and Elderfield (2007)'
Yu_Elder_df['salinity']=35


#Read Rae_2011 data
file_path=data_path/"Rae2011_calibration.csv"
Rae_df=pd.read_csv(file_path)
Rae_df['reference']='Rae et al. (2011)'

#remove all that aren't mundulus or wuellertorfi
Rae_df=Rae_df.loc[(Rae_df['species']=='Cibicidoides mundulus') | (Rae_df['species']=='Cibicidoides wuellerstorfi')]

#concatenate dataframes
Yu_Elder_Rae_df=pd.concat([Yu_Elder_df, Rae_df[['region', 'core', 'water_depth_m', 'temp_c', 
                                                'pH foram', 'DCO3_c', 'species', 'B_Ca_umolmol', 'reference', 
                                                'salinity']]], 
                          axis=0, ignore_index=True)



#convert depth to pressure
Yu_Elder_Rae_df['pressure_bar']=depth_m_to_pressure_bar(Yu_Elder_Rae_df['water_depth_m'].values)



#set Ca and Mg concentrations and salinity
calcium_conc=0.01
mag_conc=0.053


#calculate KspC and KspA from kgen
dict_ksp=kgen.calc_Ks(["KspC","KspA"], temp_c=Yu_Elder_Rae_df['temp_c'].values, sal=Yu_Elder_Rae_df['salinity'].values, p_bar=Yu_Elder_Rae_df['pressure_bar'].values, 
                       calcium=calcium_conc, magnesium=mag_conc)

Yu_Elder_Rae_df['KspC']=dict_ksp['KspC']
Yu_Elder_Rae_df['KspA']=dict_ksp['KspA']

#calculate CO3sat_c and CO3sat_c
Yu_Elder_Rae_df['CO3sat_c']=Yu_Elder_Rae_df['KspC']/calcium_conc*10**6
Yu_Elder_Rae_df['CO3sat_a']=Yu_Elder_Rae_df['KspA']/calcium_conc*10**6

#calculate omega
Yu_Elder_Rae_df['omega_c']=(Yu_Elder_Rae_df['DCO3_c']+Yu_Elder_Rae_df['CO3sat_c'])/Yu_Elder_Rae_df['CO3sat_c']
Yu_Elder_Rae_df['omega_a']=(Yu_Elder_Rae_df['DCO3_a']+Yu_Elder_Rae_df['CO3sat_a'])/Yu_Elder_Rae_df['CO3sat_a']




# import Dai data
file_path=data_path/"Brown_omega_Ca.csv"
Brown_df=pd.read_csv(file_path)
Brown_df['species']='Nuttallides umbonifera'
Brown_df['reference']='Brown et al. (2011)/Dai et al. (2023)'
Brown_df['log(omega_c)']=np.log(Brown_df['omega_c'])  


#concatenate dataframes
Yu_Elder_Rae_Brown_df=pd.concat([Yu_Elder_Rae_df, Brown_df], axis=0, ignore_index=True)


#plot B/Ca vs omega_c
fig, ax = plt.subplots()
sns.scatterplot(data=Yu_Elder_Rae_Brown_df, x='B_Ca_umolmol', y='omega_c', hue='species')
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
ax.legend(fontsize=10)
Yu_Elder_Rae_Brown_df['log(omega_c)']=np.log(Yu_Elder_Rae_Brown_df['omega_c'])


#remove nans
Yu_Elder_Rae_Brown_df.dropna(subset=['B_Ca_umolmol', 'omega_c'], inplace=True)
Yu_Elder_Rae_Brown_df.reset_index(drop=True, inplace=True)

#plot side by side with log transformed omega_c
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
sns.scatterplot(ax=ax[0], data=Yu_Elder_Rae_Brown_df, x='B_Ca_umolmol', y='omega_c', hue='species')
ax[0].set_title('raw')
ax[0].set_xlabel(r'B/Ca ($\mu$mol/mol)')
sns.scatterplot(ax=ax[1], data=Yu_Elder_Rae_Brown_df, x='B_Ca_umolmol', y='log(omega_c)', hue='species')
ax[1].set_title('transformed')
ax[1].set_xlabel(r'B/Ca ($\mu$mol/mol)')


#get species names
species_names=pd.unique(Yu_Elder_Rae_Brown_df['species'])

#make dictionaries to store linear and log fits
species_fit_dict={}
species_fit_DCO3_dict={}
xpredict_dict={}
#Make colormap for species
colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][0:len(species_names)]
species_col=dict(zip(species_names, colors))

#cycle through species and fit linear and log fits, and plot
fig, ax =plt.subplots(figsize=(6, 6))
for species in species_names:
    df=Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species]
    
    X=df['B_Ca_umolmol']
    Y=df['omega_c']
    
    xpredict=np.linspace(min(X), max(X), 30)
    xpredict_dict[species]=xpredict
    
    
    linmdl= linregress(df)
    logmdl = linregress(df, y_transform=np.log)
    
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'omega_c_linpred']=linpred(X, linmdl)
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'omega_c_logpred']=linpred(X, logmdl, y_transform=np.exp)
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'lin_r_sq']=linmdl.rsquared
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'log_r_sq']=logmdl.rsquared
    
    species_fit_dict[species]={'lin':linmdl, 'log':logmdl}
    
    species_short=species.split(' ')[1]
    
    ax.scatter(X, Y,color=species_col[species], label=species)
    ax.plot(xpredict, linpred(xpredict, linmdl), color=species_col[species], 
               label=f'r_sq = {linmdl.rsquared:.2f}')
    ax.plot(xpredict, linpred(xpredict, logmdl, y_transform=np.exp), color=species_col[species], 
               label=f'r_sq = {logmdl.rsquared:.2f}', ls='dashed')
    
    
    
    linmdl= linregress(df, target='DCO3_c')
    logmdl = linregress(df, y_transform=np.log, target='DCO3_c')
    
    species_fit_DCO3_dict={'lin':linmdl, 'log':logmdl}
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'DCO3_c_linpred']=linpred(X, linmdl)
    Yu_Elder_Rae_Brown_df.loc[Yu_Elder_Rae_Brown_df['species']==species, 'lin_DCO3_r_sq']=linmdl.rsquared
    
    
ax.legend(fontsize=8)
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
ax.set_ylabel(r'$\Omega_c$')
fig.savefig(fig_path/'B_omega_calibration.png', dpi=300)



# log fit only
fig, ax = plt.subplots()
p1=sns.scatterplot(data=Yu_Elder_Rae_Brown_df, x='B_Ca_umolmol', y='omega_c', 
                hue='species', style='reference')
for species in species_names:
    rsq=species_fit_dict[species]['log'].rsquared
    ax.plot(xpredict_dict[species], linpred(xpredict_dict[species], species_fit_dict[species]['log'], y_transform=np.exp),
            color=species_col[species], label='R$^2$ = {:.2f}'.format(rsq))
h,l = ax.get_legend_handles_labels()

for i in [4,0]:
    h.pop(i)
    l.pop(i)

p1.legend(h, l, fontsize=6, loc='upper left')
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
ax.set_ylabel(r'$\Omega_c$')
fig.savefig(fig_path/'B_omega_calibration_log.png', dpi=300)



#export the data
Yu_Elder_Rae_Brown_df.to_csv(data_path/"BCa_omega_calibration_df.csv")



## Foram data   
#import foram data
foram_path = data_path/"foram_database_240923-1419.csv"

# Read all sheets from the Excel file into a dictionary of DataFrames
foram_df = pd.read_csv(foram_path, index_col=0)


#if species_simple is wuellerstorfi, mundulus or truempyi then remove B11_corrected
foram_df.loc[foram_df['species_simple'].str.contains('wuellerstorfi|mundulus|truempyi|umbonatus'), 'B11 corrected to C. mundulus']=np.nan


#make a column designating the species to calibrate to
foram_df['B_omega_calibration']=foram_df['species_simple']

# calibrate truempyi to umbonifera
foram_df.loc[foram_df['B_omega_calibration']=='Nuttallides truempyi', 'B_omega_calibration']='Nuttallides umbonifera'


#where there is  B11_corrected, set B_omega_calibration to mundulus
foram_df.loc[pd.notnull(foram_df['B11 corrected to C. mundulus']), 'B_omega_calibration']='Cibicidoides mundulus'
#Create B11_for_omega column that uses B11 corrected to C. mundulus if available
foram_df['B11_for_omega']=foram_df['B11']
foram_df.loc[pd.notnull(foram_df['B11 corrected to C. mundulus']), 'B11_for_omega']=foram_df['B11 corrected to C. mundulus']



foram_df['omega_linfit']=np.nan
foram_df['omega_logfit']=np.nan


#cycle through species and perform linear and log predictions

for species in species_names:
    
    #perform linear and log predictions
    foram_df['omega_linfit']=np.where(foram_df['B_omega_calibration']==species, 
                                      linpred(foram_df['B11_for_omega'].values, species_fit_dict[species]['lin']), 
                                      foram_df['omega_linfit'])
    
    foram_df['omega_logfit']=np.where(foram_df['B_omega_calibration']==species,
                                        linpred(foram_df['B11_for_omega'].values, species_fit_dict[species]['log'], y_transform=np.exp),
                                        foram_df['omega_logfit'])
    



plt_df=foram_df.dropna(subset=['omega_logfit'], axis=0)
#plot
fig, ax =plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex=True)
p1=sns.scatterplot(ax=ax[0], data=plt_df, x='age_Ma', y='omega_logfit', 
                   hue='core', style='species_simple')
p2=sns.scatterplot(ax=ax[1], data=plt_df, x='age_Ma', y='B11', 
                   hue='core', style='species_simple', legend=False)
p1.legend(fontsize=6, ncols=2)

#plot
fig, ax =plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex=True)
p1=sns.scatterplot(ax=ax[0], data=plt_df, x='age_Ma', y='omega_logfit', 
                   hue='species_simple')
p2=sns.scatterplot(ax=ax[1], data=plt_df, x='age_Ma', y='B11', 
                   hue='species_simple', legend=False)
p1.legend(fontsize=6, ncols=2)





#do depth correction
foram_df['omega_linfit_0']=omega_depth_correction(foram_df['omega_linfit'], foram_df['palaeo_depth_m'])
foram_df['omega_logfit_0']=omega_depth_correction(foram_df['omega_logfit'], foram_df['palaeo_depth_m'])
foram_df['omega_linfit_3000']=omega_depth_correction(foram_df['omega_linfit'], foram_df['palaeo_depth_m'], newdepth=3000)
foram_df['omega_logfit_3000']=omega_depth_correction(foram_df['omega_logfit'], foram_df['palaeo_depth_m'], newdepth=3000)





## Temperature data


meckler_df=pd.read_csv(data_path/"Meckler2022_temp.csv")
# Create a loess fit
lowess = sm.nonparametric.lowess(meckler_df['temp_c'], meckler_df['age_Ma'], frac=0.2, it=3, delta=0.0, is_sorted=False, missing='drop', return_sorted=True)
# Create a new DataFrame with the results
lowess_df = pd.DataFrame(lowess, columns=['age_Ma', 'temp_c_smooth'])
# Merge the original data with the lowess smoothed data
merged_df = pd.merge(meckler_df, lowess_df, on='age_Ma')
# Plot the original data

fig, ax =plt.subplots(figsize=(10, 8))
ax.errorbar(merged_df['age_Ma'], merged_df['temp_c'], yerr=merged_df['err'], linestyle='None', marker='o')
# Plot the smoothed data
ax.plot(merged_df['age_Ma'], merged_df['temp_c_smooth'], color='red')
ax.set_xlabel('Age (Ma)')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.invert_xaxis()
plt.legend(['Smoothed', 'Original'])
plt.title('Meckler et al. (2022) clumped isotopes')


#interpolate onto foram data
foram_df['temp_c_Meckler']=np.nan
foram_df['temp_c_Meckler']=np.interp(foram_df['age_Ma'], meckler_df['age_Ma'], lowess_df['temp_c_smooth'])





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
ax.invert_xaxis()
plt.legend(['Calcium', 'Magnesium'])




#calculate KspC from kgen
dict_ksp=kgen.calc_Ks(["KspC"], temp_c=foram_df['temp_c_Meckler'].values, sal=35, p_bar=foram_df['pressure_bar'].values, 
                       calcium=foram_df['Ca_sw'].values/1000, magnesium=foram_df['Mg_sw'].values/1000)

foram_df['ksp_c']=dict_ksp['KspC']

#calculate [CO3]
foram_df['CO3']=foram_df['omega_logfit']*foram_df['ksp_c']/(foram_df['Ca_sw']/1000)*10**6


#plot
fig, ax =plt.subplots(figsize=(10, 6))
sns.scatterplot(data=foram_df, x='age_Ma', y='CO3', hue='core')
ax.invert_xaxis()





#isolate Cibicidoides for plotting
plt_df=foram_df.dropna(subset=['omega_logfit'], axis=0)

#plot
fig, ax =plt.subplots(nrows=2, ncols=1,  figsize=(12, 8), sharex=True, sharey=True)
p1=sns.scatterplot(ax=ax[0], data=plt_df, x='age_Ma', y='omega_logfit', hue='core', 
                style='species_simple')
p2=sns.scatterplot(ax=ax[1], data=plt_df, x='age_Ma', y='omega_logfit_3000', hue='core', 
                style='species_simple', legend=False)
p1.legend(fontsize=6, ncol=2)

#plot
fig, ax =plt.subplots(nrows=2, ncols=1,  figsize=(12, 8), sharex=True)
p1=sns.scatterplot(ax=ax[0], data=plt_df, x='age_Ma', y='B11', hue='core', 
                style='species_simple', legend=False)
p2=sns.scatterplot(ax=ax[1], data=plt_df, x='age_Ma', y='omega_logfit_3000', hue='core', 
                style='species_simple')
p2.legend(fontsize=6, ncol=2)
fig.savefig(fig_path/'B_omega3000_logfit.png', dpi=300)


#save data
timestr = time.strftime("%y%m%d-%H%M")
foram_df.to_csv(data_path/f"foram_database_BCa_{timestr}.csv")

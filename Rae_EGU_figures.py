from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import gridspec


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

#define paths
data_path=Path(os.getcwd())/"data"
fig_path=Path(os.getcwd())/"figures"


#import foram data
foram_df=pd.read_csv(data_path/"foram_dataframe.csv", index_col=0)
#import d13C d18O stacks
d13Cd18O_stack=pd.read_csv(data_path/"d13Cd18Ostack.csv")
#import d13C smooths
d13C_smooths=pd.read_csv(data_path/"d13C_smooths.csv")
#import Westerhold loess
Wester_loess_df=pd.read_csv(data_path/"Westerhold_d18Od13C_loess.csv")
#import Westerhold data
Westerhold_df=pd.read_excel(data_path/"Westerhold_full.xlsx")
#import Meckler data
Meckler_df=pd.read_csv(data_path/"Meckler2022_temp.csv")
#import BCa_calibration
BCa_cal=pd.read_csv(data_path/"BCa_omega_calibration_df.csv")
#import CO2 data
CO2_df=pd.read_csv(data_path/"Rae2022_CO2.csv", encoding='unicode_escape')
CO2_df['age']=CO2_df['age']/1000

#import Rae 2022 data
Rae2022_dict = pd.read_excel(data_path/'ea49_rae_suppl_data3.xlsx', sheet_name=None)

#import CO2 smooth
CO2_smooth=pd.read_csv(data_path/"Rae_2021_CO2_wSmooth.csv", encoding='unicode_escape')
CO2_smooth['age']=CO2_smooth['age']/1000

#import Zeebe and Tyrrell data
ZandT=pd.read_excel(data_path/"ZeebeTyrrell_2019.xlsx", sheet_name=None)


#make Cenozoic epoch boundaries dict
epoch_boundaries={'Paleocene':(65.5, 55.8), 
                  'Eocene':(55.8, 33.9), 
                  'Oligocene':(33.9, 23.0), 
                  'Miocene':(23.0, 5.3), 
                  'Pliocene':(5.3, 1.8), 
                  'Pleistocene':(1.8,0.01)}


## Figure 1: B/Ca calibration
#Make colormap for species
species_names=pd.unique(BCa_cal['species'])
colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][0:len(species_names)]


colours=['#8c2d04', '#cc4c02','#74a9cf']
custom_palette=sns.set_palette(sns.color_palette(colours))

hue_order=['Cibicidoides wuellerstorfi', 'Cibicidoides mundulus', 'Nuttallides umbonifera']
species_col=dict(zip(hue_order, colours))

#flip x and y
# plot data
fig, ax = plt.subplots()
p1=sns.scatterplot(data=BCa_cal, x='B_Ca_umolmol', y='omega_c', 
                hue='species', style='reference', 
                hue_order=hue_order, palette=custom_palette)
#plot log fits
for species in species_names:
    temp_df=BCa_cal[BCa_cal['species']==species]
    temp_df=temp_df.sort_values(by='B_Ca_umolmol')
    rsq=temp_df['log_r_sq'].values[0]
    ax.plot(temp_df['B_Ca_umolmol'], temp_df['omega_c_logpred'], 
            color=species_col[species], label='R$^2$ = {:.2f}'.format(rsq))
#remove unwanted legend entries
h,l = ax.get_legend_handles_labels()
for i in [4,0]:
    h.pop(i)
    l.pop(i)
p1.legend(h, l, fontsize=6, loc='upper left')
ax.set_xlabel(r'B/Ca ($\mu$mol/mol)')
ax.set_ylabel(r'$\Omega_c$')
fig.savefig(fig_path/'RaeEGU_B_omega_calibration_log.png', dpi=300)




## Figure 2: CO2
fig, ax = plt.subplots(figsize=(6, 6))
p1=sns.scatterplot(data=CO2_df, x='age', y='xco2')
ax.invert_xaxis()
ax.set_xlabel('Age (Ma)')
ax.set_ylabel('CO$_2$ (ppm)')
fig.savefig(fig_path/'RaeEGU_CO2.png', dpi=300)






## Figure 3: d13C d18O stacks
d13Cd18O_stack['ocean_basin']='Other'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('Atl'),'ocean_basin']='Atlantic'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('Pac'),'ocean_basin']='Pacific/Indian'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('SOcn'),'ocean_basin']='Southern'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('IOcn'),'ocean_basin']='Pacific/Indian'


from loess.loess_1d import loess_1d

for ocean in ['Atlantic', 'Pacific/Indian', 'Southern']:
    temp_df=d13Cd18O_stack[d13Cd18O_stack['ocean_basin']==ocean]

    xout, yout, wout = loess_1d(temp_df['Age'].values, temp_df['d13C'].values, 
                            xnew=None, degree=1, npoints=100, rotate=False, sigy=None)

    
    d13Cd18O_stack.loc[d13Cd18O_stack['ocean_basin']==ocean, 'd13C_smooth']=yout
    


d13C_smooths_melt=pd.melt(d13C_smooths, id_vars=['age_Ma'], 
                             value_vars=['Atlantic', 'Pacific', 'Southern'], 
                             var_name='ocean_basin', value_name='d13C_smooth')

d13Cd18O_stack.sort_values(by='Age', inplace=True)

d18Oloess=d13Cd18O_stack[['Age', 'd18O']].copy()
#drop na
d18Oloess.dropna(inplace=True)
#too many Holocene data points, remove some
d18Oloess=d18Oloess.iloc[170:]

d18Ox, d18Oy, wout = loess_1d(d18Oloess['Age'].values, d18Oloess['d18O'].values, 
                            xnew=np.arange(72), degree=1, npoints=500, rotate=False, sigy=None)



#plot d13C and d18O with 1Ma smoothing
fig, ax = plt.subplots(nrows=2, figsize=(10, 6), sharex=True)
p1=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d18O', hue='ocean_basin', ax=ax[0], 
                   edgecolor='none', s=1, palette='pastel', alpha=0.3)
p1b=ax[0].plot(d18Ox, d18Oy, color='black', label='All data smooth')
p2=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d13C', hue='ocean_basin', ax=ax[1], 
                   edgecolor='none', s=1, palette='pastel', alpha=0.3)
p3=sns.lineplot(data=d13C_smooths_melt, x='age_Ma', y='d13C_smooth', hue='ocean_basin', 
                ax=ax[1], palette='dark', hue_order=['Atlantic', 'Pacific', 'Southern'], legend=False)
ax[0].legend(markerscale=5, fontsize=6, loc='lower left')
ax[1].legend(markerscale=5, fontsize=6, loc='lower left')
ax[0].set_ylim(-2, 6)
ax[1].set_ylim(-2, 3)
ax[0].set_xlim(0, 65)
ax[0].invert_xaxis()
ax[0].invert_yaxis()


#plot d13C and d18O with 1Ma smoothing
fig, ax = plt.subplots(nrows=2, figsize=(6, 6), sharex=True)
p1=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='black', alpha=0.3)
p1b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='red', label='Westerhold et al. (2021) loess', 
                 ax=ax[0])
p2=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d13C', hue='ocean_basin', ax=ax[1], 
                   edgecolor='none', s=1, palette='pastel', alpha=0.3)
p3=sns.lineplot(data=d13C_smooths_melt, x='age_Ma', y='d13C_smooth', hue='ocean_basin', 
                ax=ax[1], palette='dark', hue_order=['Atlantic', 'Pacific', 'Southern'], legend=False)
ax[0].legend(markerscale=5, fontsize=8, loc='upper right')
ax[1].legend(markerscale=5, fontsize=8, loc='lower left')
ax[0].set_ylim(-2, 6)
ax[1].set_ylim(-2, 3)
ax[0].set_xlim(0, 65)
ax[0].invert_xaxis()
ax[0].invert_yaxis()
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.3)
#remove the box around the plots

ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['left'].set_visible(False)
#put the y axis labels on the right for the bottom plot
ax[1].yaxis.tick_right()
#and the label on the right
ax[1].yaxis.set_label_position('right')
#remove the background colour
ax[0].set_facecolor('none')
ax[1].set_facecolor('none')











"""
#plot d13C and d18O with 100 data point smoothing
fig, ax = plt.subplots(nrows=2, figsize=(10, 6), sharex=True)
p1=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d18O', hue='ocean_basin', ax=ax[0], 
                   edgecolor='none', s=1, palette='pastel', alpha=0.3)
p2=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d13C', hue='ocean_basin', ax=ax[1], 
                   edgecolor='none', s=1, palette='pastel', alpha=0.3)
p3=sns.lineplot(data=d13Cd18O_stack, x='Age', y='d13C_smooth', hue='ocean_basin', 
                ax=ax[1], palette='dark')
ax[0].legend(markerscale=5, fontsize=6, loc='upper right')
ax[1].legend(markerscale=5, fontsize=6, loc='lower right')
ax[0].set_ylim(-2, 6)
ax[1].set_ylim(-2, 3)
ax[0].invert_xaxis()
ax[0].invert_yaxis()
"""





## Figure 4: Temperature


from loess.loess_1d import loess_1d
temp_df=pd.read_csv(data_path/"Meckler2022_temp.csv")

xout, yout, wout = loess_1d(temp_df['age_Ma'].values, temp_df['temp_c'].values, 
                            xnew=None, degree=1, npoints=8, rotate=False, sigy=None)

#plot Meckler data
fig, ax =plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=temp_df, x='age_Ma', y='temp_c', color='black', label='Meckler et al. (2022) data')
p2=ax.plot(xout, yout, color='red', label='Loess fit (n=8)')
ax.set_xlim(0, 65)
ax.invert_xaxis()
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlabel('Age (Ma)')
ax.legend()



## B/Ca

#subset foram_df to only include data with omega_c
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']

custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order)
ax.legend(fontsize=8, ncol=2)
ax.set(ylabel='B/Ca ($\mu$mol/mol)')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()





BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

foram_df_B11_addedtime0=pd.concat([foram_df_B11, BCa_cal], ignore_index=True, axis=0)

fig, ax = plt.subplots(figsize=(6, 3))
p1=sns.scatterplot(data=foram_df_B11_addedtime0, x='age_Ma', y='B11', 
                hue='ocean_basin', palette='Set2', s=10)
ax.legend(fontsize=8, ncol=2)
ax.set(ylabel='B/Ca ($\mu$mol/mol)')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()


## Figure 5: Omega


fig, ax = plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order)
ax.legend(fontsize=8, ncol=2, loc='upper left')
ax.set(ylabel=r'$\Omega_c$')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()



fig, ax = plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette='Set2')
ax.legend(fontsize=8, ncol=2, loc='upper left')
ax.set(ylabel=r'$\Omega_c$')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()



## Figure 6: [CO3-2]

fig, ax = plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order)
ax.legend(fontsize=8, ncol=2, loc='upper left')
ax.set(ylabel='[CO$_3^{2-}$] ($\mu$mol/kg)')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()




fig, ax = plt.subplots(figsize=(12, 6))
p1=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='ocean_basin', palette='Set2')
ax.legend(fontsize=8, ncol=2, loc='upper left')
ax.set(ylabel='[CO$_3^{2-}$] ($\mu$mol/kg)')
ax.set(xlabel='Age (Ma)')
ax.invert_xaxis()








## Fig 1: Westerhold, pH, CO2

fig, ax = plt.subplots(figsize=(5, 7), nrows=3, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False)
p2=sns.scatterplot(data=CO2_df, x='age', y='pH', ax=ax[1], color='darkorchid', 
                   label='Boron-based pH', s=15, legend=False)
p3_a=sns.scatterplot(data=CO2_df, x='age', y='xco2', ax=ax[2], color='orange', 
                   label='Boron-based CO$_2$', s=25, legend=False, marker='^')
p3_b=sns.lineplot(data=CO2_smooth, x='age', y='xco2_smooth', ax=ax[2], color='darkgoldenrod', 
                  label='Smoothed CO$_2$', legend=False, linestyle='--')

ax[0].set_ylim(-1.5, 5.5)
#set y ticks
ax[0].set_yticks(np.arange(0, 6, 2))
ax[2].set_yticks(np.arange(0, 3000, 500))
ax[0].set_xlim(0, 70)

ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[1].set_ylabel(r'Surface ocean pH')
ax[2].set_ylabel(r'Atmospheric CO$_2$ (ppm)')
ax[2].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',)
#remove x ticks on all but the bottom plot
fig.savefig(fig_path/'EGU_fig1.png', dpi=300)




## Fig 2: Westerhold, CCD, Omega, Ca, CO3 (surface), CO2

#import CCD data
Lyle_CCD=pd.read_excel(data_path/"CCD"/"Lyle_2008_CCD_BoudreauGC 1.xlsx")
Tyrell_CCD=pd.read_excel(data_path/"CCD"/"Tyrrell_2004_CCD_GC 1.xlsx")
VanAndel_CCD=pd.read_excel(data_path/"CCD"/"VanAndel_1975_CCD_TyrrellGC 1.xlsx")
Palike2012=pd.read_csv(data_path/"CCD"/"Palike2012.csv")
Palike2012['CCDdepth']=Palike2012['DepthEqCCD']/1000
Palike2012['age']=Palike2012['age']/1000
Palike2012['reference']='Palike'
VanAndel_CCD['reference']='Van Andel'
Tyrell_CCD['reference']='Tyrrell'
Lyle_CCD['reference']='Lyle'
VanAndel_CCD.sort_values(by='age', inplace=True)
CCD_df=pd.concat([Lyle_CCD, Tyrell_CCD, VanAndel_CCD, Palike2012[['age', 'CCDdepth', 'reference']]], ignore_index=True, axis=0)
CCD_df.dropna(inplace=True,axis=0)

#Import old Ca data
Ca_old=pd.read_csv(data_path/"Ca_old.csv")
Ca_old['err_up']=Ca_old['HalCaUP']+Ca_old['HalCa']
Ca_old['err_low']=Ca_old['HalCa']-Ca_old['HalCaLOW']



fig, ax = plt.subplots(figsize=(6, 10), nrows=6, sharex=True, constrained_layout=False)

#d18O
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False)
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))
ax[0].set_ylabel(r'Benthic $\delta^{18}$O')
ax[0].invert_yaxis()

#CCD
p2=sns.lineplot(data=CCD_df, x='age', y='CCDdepth', ax=ax[1], style='reference', )
ax[1].set_ylabel('CCD depth (km)')
ax[1].invert_yaxis()
ax[1].legend(loc='center left', bbox_to_anchor=(-0.2, 0.5))

#Omega (surface)
p3=[]
colours=['#1b9e77', '#d95f02', '#7570b3', '#e7298a']
df_names=['CO2system_pH_omega65', 'CO2system_pH_alkalinity', 'CO2system_pH_dic', 'CO2system_pH_ccd']
df_labels=[r'$\Omega = 6.5$', 'Alkalinity-based', 'DIC-based', 'CCD-based']
for i, (p3_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p3_df].copy()
    temp_df['age']=temp_df['age']/1000
    p3.append(sns.lineplot(data=temp_df, x='age', y='saturation_state',ax=ax[2], 
                       color=c, label=lab, legend=False))
ax[2].set_ylim(2, 14)
ax[2].set_ylabel('Surface $\Omega$')


#Ca
p4=ax[3].errorbar(x=Ca_old['HalAbsAge'], y=Ca_old['HalCa'], 
                  yerr=[Ca_old['HalCaLOW'], Ca_old['HalCaUP']], fmt='D', 
                  color='grey', markersize=5, label='Halite Ca')
ax[3].set_ylim(0, 25)
ax[3].set_ylabel(r'[Ca$^{2+}$] (mmol/kg)')

#CO3
p5=[]
for i, (p5_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p5_df].copy()
    temp_df['age']=temp_df['age']/1000
    p5.append(sns.lineplot(data=temp_df, x='age', y='co3',ax=ax[4], 
                       color=c, label=lab))
ax[4].set_ylabel(r'[CO$_3^{2-}$] ($\mu$mol/kg)')
ax[4].legend(loc='upper left', bbox_to_anchor=(0, 1.2))

#CO2
p6_a=sns.scatterplot(data=CO2_df, x='age', y='xco2', ax=ax[5], color='orange', 
                label='Boron-based CO$_2$', s=25, legend=False, marker='^')
p6_b=sns.lineplot(data=CO2_smooth, x='age', y='xco2_smooth', ax=ax[5], color='darkgoldenrod', 
                  label='Smoothed CO$_2$', legend=False, linestyle='--')
ax[5].set_ylabel(r'Atmospheric CO$_2$ (ppm)')



ax[0].set_xlim(0, 70)
ax[5].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()



#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.2)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',)
#remove x ticks on all but the bottom plot

fig.savefig(fig_path/'EGU_fig2.png', dpi=300)





## Fig 3: B/Ca calibration


#Make colormap for species
species_names=pd.unique(BCa_cal['species'])
colours=['#8c2d04', '#cc4c02','#74a9cf']
custom_palette=sns.set_palette(sns.color_palette(colours))

hue_order=['Cibicidoides wuellerstorfi', 'Cibicidoides mundulus', 'Nuttallides umbonifera']
species_col=dict(zip(hue_order, colours))


# plot data
fig, ax = plt.subplots(figsize=(12, 5), ncols=2, sharey=True)
p1=sns.scatterplot(data=BCa_cal, y='B_Ca_umolmol', x='omega_c', 
                hue='species', style='reference', 
                hue_order=hue_order, palette=custom_palette, ax=ax[0])
#plot log fits
for species in species_names:
    temp_df=BCa_cal[BCa_cal['species']==species]
    temp_df=temp_df.sort_values(by='B_Ca_umolmol')
    rsq=temp_df['log_r_sq'].values[0]
    ax[0].plot(temp_df['omega_c_logpred'], temp_df['B_Ca_umolmol'],
            color=species_col[species], label='R$^2$ = {:.2f}'.format(rsq))
#remove unwanted legend entries
h,l = ax[0].get_legend_handles_labels()
for i in [4,0]:
    h.pop(i)
    l.pop(i)
ax[0].legend(h, l, fontsize=6, loc='lower right')
ax[0].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[0].set_xlabel(r'$\Omega_c$')


p2=sns.scatterplot(data=BCa_cal, y='B_Ca_umolmol', x='DCO3_c', 
                hue='species', style='reference', 
                hue_order=hue_order, palette=custom_palette, ax=ax[1])
#plot log fits
for species in species_names:
    temp_df=BCa_cal[BCa_cal['species']==species]
    temp_df=temp_df.sort_values(by='B_Ca_umolmol')
    rsq=temp_df['lin_DCO3_r_sq'].values[0]
    ax[1].plot(temp_df['DCO3_c_linpred'], temp_df['B_Ca_umolmol'],
            color=species_col[species], label='R$^2$ = {:.2f}'.format(rsq))
#remove unwanted legend entries
h,l = ax[1].get_legend_handles_labels()
for i in [4,0]:
    h.pop(i)
    l.pop(i)
ax[1].legend(h, l, fontsize=6, loc='lower right')
ax[1].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1].set_xlabel(r'$\Delta$[CO$_3^{2-}$] ($\mu$mol/kg)')
plt.subplots_adjust(wspace=0)



fig.savefig(fig_path/'EGU_fig3.png', dpi=300)





## Fig 4a: Westerhold and B/Ca

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']

custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(6, 6), nrows=2, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, 
                  zorder=3)
ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1].set_ylim(25, 250)
ax[1].legend(fontsize=8, ncol=2, loc='lower left')

ax[0].set_xlim(0, 70)
ax[1].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)



fig.savefig(fig_path/'EGU_fig4a.png', dpi=300)




## Fig 4b_box: Westerhold and B/Ca w modern

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]



BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']
BCa_cal['species_simple']=BCa_cal['species']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

foram_df_B11_addedtime0=pd.concat([foram_df_B11, BCa_cal], ignore_index=True, axis=0)



colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']


custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, ncol=2, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



custom_palette_small=sns.set_palette(sns.color_palette(['#74a9cf', '#cc4c02', '#8c2d04']))

p3=sns.boxplot(data=BCa_cal, x='species_simple', hue='species_simple', y='B11', ax=ax[1,1], 
               palette=custom_palette_small, legend=False, whis=[0, 100], 
               hue_order=['Nuttallides umbonifera', 'Cibicidoides mundulus', 'Cibicidoides wuellerstorfi'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')

fig.savefig(fig_path/'EGU_fig4b_box.png', dpi=300)






## Fig 4b_strip: Westerhold and B/Ca w modern

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]



BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']
BCa_cal['species_simple']=BCa_cal['species']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

foram_df_B11_addedtime0=pd.concat([foram_df_B11, BCa_cal], ignore_index=True, axis=0)



colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']


custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, ncol=2, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



custom_palette_small=sns.set_palette(sns.color_palette(['#74a9cf', '#cc4c02', '#8c2d04']))


p3b=sns.stripplot(data=BCa_cal, x='species_simple',  hue='species_simple', y='B11', ax=ax[1,1], 
               palette=custom_palette_small, legend=False, jitter=0.4, dodge=True, 
               size=3, alpha=0.8, 
               hue_order=['Nuttallides umbonifera', 'Cibicidoides mundulus', 'Cibicidoides wuellerstorfi'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')

fig.savefig(fig_path/'EGU_fig4b_strip.png', dpi=300)











## Fig 5a: Westerhold and B/Ca (basin)

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(6, 6), nrows=2, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, 
                  zorder=3)
ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='ocean_basin', palette=custom_palette,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1].set_ylim(25, 225)
ax[1].legend(fontsize=10, loc='lower left')

ax[0].set_xlim(0, 70)
ax[1].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)



fig.savefig(fig_path/'EGU_fig5a.png', dpi=300)




## Fig 5b_box: Westerhold and B/Ca (basin) modern box






colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='ocean_basin', palette=custom_palette, 
                ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



p3b=sns.boxplot(data=BCa_cal, x='ocean_basin',  hue='ocean_basin', y='B11', ax=ax[1,1], 
               palette=custom_palette, legend=False, whis=[0, 100], 
               order=['Atlantic', 'Southern', 'Pacific'],
               hue_order=['Atlantic', 'Southern', 'Pacific'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')


fig.savefig(fig_path/'EGU_fig5b_box.png', dpi=300)

## Fig 5b_strip: Westerhold and B/Ca (basin) modern box strip

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11', 
                hue='ocean_basin', palette=custom_palette, 
                ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



p3b=sns.stripplot(data=BCa_cal, x='ocean_basin',  hue='ocean_basin', y='B11', ax=ax[1,1], 
               palette=custom_palette, legend=False, jitter=0.4, dodge=True, 
               size=3, alpha=0.8, 
               order=['Atlantic', 'Southern', 'Pacific'],
               hue_order=['Atlantic', 'Southern', 'Pacific'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')


fig.savefig(fig_path/'EGU_fig5b_strip.png', dpi=300)






## Fig 6a: Westerhold and species-corrected B/Ca

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']

custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(6, 6), nrows=2, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, 
                  zorder=3)
ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1].set_ylim(25, 250)
ax[1].legend(fontsize=8, ncol=2, loc='lower left')

ax[0].set_xlim(0, 70)
ax[1].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)



fig.savefig(fig_path/'EGU_fig6a.png', dpi=300)




## Fig 6b_box: Westerhold and B/Ca w modern

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]



BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']
BCa_cal['species_simple']=BCa_cal['species']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

foram_df_B11_addedtime0=pd.concat([foram_df_B11, BCa_cal], ignore_index=True, axis=0)



colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']


custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, ncol=2, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



custom_palette_small=sns.set_palette(sns.color_palette(['#74a9cf', '#cc4c02', '#8c2d04']))

p3=sns.boxplot(data=BCa_cal, x='species_simple', hue='species_simple', y='B11', ax=ax[1,1], 
               palette=custom_palette_small, legend=False, whis=[0, 100], 
               hue_order=['Nuttallides umbonifera', 'Cibicidoides mundulus', 'Cibicidoides wuellerstorfi'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')

fig.savefig(fig_path/'EGU_fig6b_box.png', dpi=300)






## Fig 6b_strip: Westerhold and B/Ca w modern

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]



BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']
BCa_cal['species_simple']=BCa_cal['species']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

foram_df_B11_addedtime0=pd.concat([foram_df_B11, BCa_cal], ignore_index=True, axis=0)



colours=['#045a8d', '#74a9cf', 
         '#ffffd4', 
         '#fee391', 
         '#fec44f',
         '#fe9929', 
         '#ec7014',
         '#cc4c02',
         '#8c2d04']
hue_order=['Nuttallides truempyi', 'Nuttallides umbonifera', 
           'Cibicidoides grimsdalei',
           'Cibicidoides havanensis', 
           'Cibicidoides subspiratus',
           'Cibicidoides eocaenus',
           'Cibicidoides laurisae', 
           'Cibicidoides mundulus',
           'Cibicidoides wuellerstorfi']


custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='species_simple', palette=custom_palette, 
                hue_order=hue_order,ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, ncol=2, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



custom_palette_small=sns.set_palette(sns.color_palette(['#74a9cf', '#cc4c02', '#8c2d04']))


p3b=sns.stripplot(data=BCa_cal, x='species_simple',  hue='species_simple', y='B11', ax=ax[1,1], 
               palette=custom_palette_small, legend=False, jitter=0.4, dodge=True, 
               size=3, alpha=0.8, 
               hue_order=['Nuttallides umbonifera', 'Cibicidoides mundulus', 'Cibicidoides wuellerstorfi'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')

fig.savefig(fig_path/'EGU_fig6b_strip.png', dpi=300)











## Fig 7a: Westerhold and B/Ca (basin)

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(6, 6), nrows=2, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, 
                  zorder=3)
ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='ocean_basin', palette=custom_palette,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1].set_ylim(25, 225)
ax[1].legend(fontsize=10, loc='lower left')

ax[0].set_xlim(0, 70)
ax[1].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)



fig.savefig(fig_path/'EGU_fig7a.png', dpi=300)




## Fig7b_box: Westerhold and B/Ca (basin) modern box






colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='ocean_basin', palette=custom_palette, 
                ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



p3b=sns.boxplot(data=BCa_cal, x='ocean_basin',  hue='ocean_basin', y='B11', ax=ax[1,1], 
               palette=custom_palette, legend=False, whis=[0, 100], 
               order=['Atlantic', 'Southern', 'Pacific'],
               hue_order=['Atlantic', 'Southern', 'Pacific'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')


fig.savefig(fig_path/'EGU_fig7b_box.png', dpi=300)

## Fig7b_strip: Westerhold and B/Ca (basin) modern box strip

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

BCa_cal['age_Ma']=0
BCa_cal['B11']=BCa_cal['B_Ca_umolmol']
BCa_cal['ocean_basin']=BCa_cal['region']

ocean_basins={'Norwegian':'Atlantic', 'Canaries':'Atlantic', 'Walvis East':'Atlantic',
              'Namibia':'Atlantic','N Atlantic':'Atlantic', 'Walvis West':'Atlantic','Indian':'Pacific', 
              'Atlantic':'Atlantic', 'Southern':'Southern', 'Pacific':'Pacific'}

BCa_cal['ocean_basin']=BCa_cal['region'].map(ocean_basins)   

fig, ax = plt.subplots(figsize=(7, 6), nrows=2, ncols=2 ,
                       constrained_layout=False, gridspec_kw={'width_ratios': [6, 1]})

p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0,0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0,0], legend=False, 
                  zorder=3)
ax[0,0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0,0].set_ylim(-1.5, 5.5)
ax[0,0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='B11_for_omega', 
                hue='ocean_basin', palette=custom_palette, 
                ax=ax[1,0],zorder=3, s=20)
ax[1,0].set_ylabel(r'Species-corrected B/Ca ($\mu$mol/mol)')
ax[1,0].set_ylim(25, 275)
ax[1,0].legend(fontsize=8, loc='lower left')

ax[0,0].set_xlim(0, 70)
ax[1,0].set_xlim(0, 70)
ax[0,0].set_xlabel('')
ax[1,0].set_xlabel('Age (Ma)')
ax[0,0].invert_xaxis()
ax[0,0].invert_yaxis()
ax[1,0].invert_yaxis()
ax[1,0].invert_xaxis()



p3b=sns.stripplot(data=BCa_cal, x='ocean_basin',  hue='ocean_basin', y='B11', ax=ax[1,1], 
               palette=custom_palette, legend=False, jitter=0.4, dodge=True, 
               size=3, alpha=0.8, 
               order=['Atlantic', 'Southern', 'Pacific'],
               hue_order=['Atlantic', 'Southern', 'Pacific'])
ax[1,1].set_ylim(25, 275)
ax[1,1].invert_yaxis()

ax[1,1].axis('off')
ax[1,1].text(1, 60, s='Modern',  size='medium', color='black', horizontalalignment='center')
#ax[1,1].title.set_text('Modern')
#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.4)
plt.subplots_adjust(wspace=0.3)
#remove the box around the plots

for i, axis in enumerate(ax[:,0]):
    #if even
    if i%2 == 0:
        axis.spines['right'].set_visible(False)
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position('right')
    if i > 0:
        axis.spines['top'].set_visible(False)
    if i < len(ax[:,0])-1:
        axis.spines['bottom'].set_visible(False)
        plt.setp(axis.get_xticklabels(), visible=False)
        plt.setp(axis.get_xticklines(), visible=False)

    axis.set_facecolor('none')
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)


ax[0, 1].axis('off')


fig.savefig(fig_path/'EGU_fig7b_strip.png', dpi=300)




## Fig 8a: Westerhold and Omega (basin)


foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

fig, ax = plt.subplots(figsize=(6, 6), nrows=2, sharex=True, constrained_layout=False)
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False,
                   zorder=2)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, 
                  zorder=3)
ax[0].set_ylabel(r'Benthic foraminiferal $\delta^{18}$O')
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette=custom_palette,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'$\Omega_c$')
#ax[1].set_ylim(25, 225)
ax[1].legend(fontsize=10, loc='center left')

ax[0].set_xlim(0, 70)
ax[1].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()
ax[0].invert_yaxis()
#ax[1].invert_yaxis()

#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.5)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--',
                     zorder=1)



fig.savefig(fig_path/'EGU_fig8.png', dpi=300)



## Fig 9: Westerhold, Omega (basin), Ca, CO3

foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]
Ca_df=pd.read_excel(data_path/"calcium_magnesium.xlsx", sheet_name='calcium')



fig, ax = plt.subplots(figsize=(6, 10), nrows=4, sharex=True, constrained_layout=False)

#d18O
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False, zorder=3)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, zorder=4)
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))
ax[0].set_ylabel(r'Benthic $\delta^{18}$O')
ax[0].set_zorder(4)
ax[0].invert_yaxis()

#Omega (benthic)
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette=custom_palette,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'$\Omega_c$')
#ax[1].set_ylim(25, 225)
ax[1].legend(fontsize=10, loc='center left')
ax[1].set_ylim(0.5, 2)


#Ca

x=Ca_df['age'].values
y=Ca_df['median'].values
upper=Ca_df['95% quantile'].values
lower=Ca_df['5% quantile'].values
p3=ax[2].plot(x, y, color='black', zorder=3)
p3err=ax[2].fill_between(x, lower, upper, color='#d08f8d', alpha=0.7, zorder=3)
ax[2].set_ylabel(r'[Ca$^{2+}$] (mmol/kg)')
ax[2].set_ylim(8, 42)

#CO3

p4=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='ocean_basin', palette=custom_palette, ax=ax[3], s=20, zorder=3)
ax[3].legend(fontsize=10, loc='center left')
ax[3].set(ylabel='[CO$_3^{2-}$] ($\mu$mol/kg)')
ax[3].set(xlabel='Age (Ma)')
ax[3].invert_xaxis()
ax[3].set_ylim(20, 100)



ax[0].set_xlim(0, 70)
ax[3].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()



#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.3)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--', zorder=1)
#remove x ticks on all but the bottom plot

fig.savefig(fig_path/'EGU_fig9.png', dpi=300)



## Fig 10: Westerhold, Omega (basin), CCD (basin), CO3 (basin), d13C (basin)

d13Cd18O_stack['ocean_basin']='Other'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('Atl'),'ocean_basin']='Atlantic'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('Pac'),'ocean_basin']='Pacific/Indian'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('SOcn'),'ocean_basin']='Southern'
d13Cd18O_stack.loc[d13Cd18O_stack['Location'].str.contains('IOcn'),'ocean_basin']='Pacific/Indian'


d13C_smooths_melt=pd.melt(d13C_smooths, id_vars=['age_Ma'], 
                             value_vars=['Atlantic', 'Pacific', 'Southern'], 
                             var_name='ocean_basin', value_name='d13C_smooth')

d13Cd18O_stack.sort_values(by='Age', inplace=True)
foram_df_B11=foram_df.loc[pd.notnull(foram_df['omega_logfit'])]

fig, ax = plt.subplots(figsize=(6, 10), nrows=5, sharex=True, constrained_layout=False)

#d18O
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False, zorder=3)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, zorder=4)
ax[0].set_ylim(-1.5, 5.5)
ax[0].set_yticks(np.arange(0, 6, 2))
ax[0].set_ylabel(r'Benthic $\delta^{18}$O')
ax[0].invert_yaxis()
ax[0].set_zorder(4)


#Omega (benthic)
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

p2=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette=custom_palette,ax=ax[1],zorder=3, s=20)
ax[1].set_ylabel(r'$\Omega_c$')
#ax[1].set_ylim(25, 225)
ax[1].legend(fontsize=10, bbox_to_anchor=(0.2, 0.9), loc='center left')
ax[1].set_ylim(0.3, 2)
ax[1].set_yticks([0.5, 1, 1.5, 2])
ax[1].set_zorder=3

#CCD
CCD_basin_dict={'Palike':'Pacific', 'Lyle':'Pacific'}

CCD_df['ocean_basin']='Global'
CCD_df.loc[CCD_df['reference'].isin(CCD_basin_dict.keys()), 'ocean_basin']=CCD_df.loc[CCD_df['reference'].isin(CCD_basin_dict.keys()), 'reference'].map(CCD_basin_dict)
colours={'Palike': '#7570b3', 'Lyle': '#7570b3', 'Tyrrell': '#2DC7FF', 'Van Andel': '#2DC7FF'}
styles={'Palike': 'dashdot', 'Lyle': 'solid', 'Tyrrell': 'dashed', 'Van Andel': 'dotted'}


for key, value in colours.items():
    ax[2].plot(CCD_df.loc[CCD_df['reference']==key, 'age'], CCD_df.loc[CCD_df['reference']==key, 'CCDdepth'], 
               color=value, linestyle=styles[key], 
               label=key+' ('+CCD_df.loc[CCD_df['reference']==key, 'ocean_basin'].values[0]+')', 
               zorder=3)
ax[2].set_zorder(3)
ax[2].set_ylabel('CCD depth (km)')
ax[2].set_ylim(2.5, 5.5)
ax[2].set_yticks([3, 4, 5])
ax[2].invert_yaxis()
ax[2].legend(loc='lower left', bbox_to_anchor=(0, -0.1))



#CO3
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p4=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='ocean_basin', palette=custom_palette, ax=ax[3], s=20, 
                legend=False, zorder=3)
#ax[3].legend(fontsize=10, loc='center left')
ax[3].set(ylabel='[CO$_3^{2-}$] ($\mu$mol/kg)')
ax[3].set(xlabel='Age (Ma)')
ax[3].set_yticks([25, 50, 75, 100])
ax[3].invert_xaxis()
ax[3].set_ylim(5, 110)


#d13C
hue_order=['Atlantic', 'Southern', 'Pacific/Indian']
colours=['#7fc97f', '#fdc086', '#beaed4']
custom_palette=sns.set_palette(sns.color_palette(colours))

p5a=sns.scatterplot(data=d13Cd18O_stack, x='Age', y='d13C', hue='ocean_basin', ax=ax[4], 
                   edgecolor='none', s=1, palette=custom_palette, alpha=0.2, 
                   hue_order=hue_order, legend=False, zorder=3)

hue_order=['Atlantic', 'Southern', 'Pacific']
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))

p5b=sns.lineplot(data=d13C_smooths_melt, x='age_Ma', y='d13C_smooth', hue='ocean_basin', 
                ax=ax[4], palette=custom_palette, hue_order=hue_order, zorder=4)
ax[4].legend(loc='lower left', ncols=3)
ax[4].set_ylim(-1.8, 3)
ax[4].set_ylabel(r'Benthic $\delta^{13}$C')


ax[0].set_xlim(0, 70)
ax[4].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()



#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.3)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--', zorder=1)
#remove x ticks on all but the bottom plot

fig.savefig(fig_path/'EGU_fig10.png', dpi=300)



## Fig 11a: Westerhold, Omega (basin), 






CCD_CO3_interp=np.interp(ZandT['boron']['Age(Ma)'].values, ZandT['CO3model']['age'].values, ZandT['CO3model']['CO3_CCD'].values)

CO3_offset=ZandT['boron']['[CO32-](mumol/kg)'].values-CCD_CO3_interp

max_age=ZandT['boron']['Age(Ma)'].max()

foram_df_B11_agecut=foram_df_B11.loc[foram_df_B11['age_Ma']<max_age]

foram_df_B11_agecut['surface_deep_CO3_offset']=np.interp(foram_df_B11_agecut['age_Ma'].values, ZandT['boron']['Age(Ma)'].values, CO3_offset)

foram_df_B11_agecut['surface_CO3']=foram_df_B11_agecut['CO3']+foram_df_B11_agecut['surface_deep_CO3_offset']









fig, ax = plt.subplots(figsize=(6, 10), nrows=5, sharex=True, constrained_layout=False)

#d18O
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False)
ax[0].set_ylim(-1.5, 6)
ax[0].set_yticks(np.arange(0, 6, 2))
ax[0].set_ylabel(r'Benthic $\delta^{18}$O')
ax[0].invert_yaxis()
ax[0].set_zorder(3)


#Omega (surface)
p2=[]
colours=['#FCB0B3', '#F93943', '#7EB2DD', '#445E93']
df_names=['CO2system_pH_omega65', 'CO2system_pH_alkalinity', 'CO2system_pH_dic', 'CO2system_pH_ccd']
df_labels=[r'$\Omega = 6.5$', 'Alkalinity-based', 'DIC-based', 'CCD-based']
for i, (p3_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p3_df].copy()
    temp_df['age']=temp_df['age']/1000
    p2.append(sns.lineplot(data=temp_df, x='age', y='saturation_state',ax=ax[1], 
                       color=c, label=lab, legend=False))
ax[1].set_ylim(3, 9)
ax[1].set_ylabel('Surface $\Omega_c$')



#Omega (benthic)
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p3=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette=custom_palette,ax=ax[2],zorder=3, s=20)
ax[2].set_ylabel(r'Deep $\Omega_c$')
#ax[1].set_ylim(25, 225)
ax[2].legend(fontsize=10, loc='upper left')
ax[2].set_ylim(0.3, 2)
ax[2].set_yticks([0.5, 1, 1.5, 2])



#CO3 (surface)
p4=[]
colours=['#FCB0B3', '#F93943', '#7EB2DD', '#445E93']
for i, (p5_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p5_df].copy()
    temp_df['age']=temp_df['age']/1000
    p4.append(sns.lineplot(data=temp_df, x='age', y='co3',ax=ax[3], 
                       color=c, label=lab))
ax[3].set_ylabel(r'Surface [CO$_3^{2-}$] ($\mu$mol/kg)')
ax[3].legend(loc='upper left')


#CO3 (benthic)
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p4=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='ocean_basin', palette=custom_palette, ax=ax[4], s=20, 
                legend=False)
#ax[3].legend(fontsize=10, loc='center left')
ax[4].set(ylabel='Deep [CO$_3^{2-}$] ($\mu$mol/kg)')
ax[4].set(xlabel='Age (Ma)')
ax[4].set_yticks([25, 50, 75, 100])
ax[4].invert_xaxis()
ax[4].set_ylim(5, 110)




ax[0].set_xlim(0, 70)
ax[4].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()



#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.2)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--', zorder=1)
#remove x ticks on all but the bottom plot

fig.savefig(fig_path/'EGU_fig11a.png', dpi=300)


## Fig 11b: Westerhold, Omega (basin), 





fig, ax = plt.subplots(figsize=(6, 10), nrows=5, sharex=True, constrained_layout=False)

#d18O
p1_a=sns.scatterplot(data=Westerhold_df, x='Age', y='d18O_adj', ax=ax[0], 
                   edgecolor='none', s=1, color='grey', alpha=0.5, 
                   label='Westerhold et al. (2020) data', legend=False, zorder=3)
p1_b=sns.lineplot(data=Wester_loess_df, x='age_Ma', y='loess_long_d18O', color='black', 
                  label='Westerhold et al. (2020) loess', ax=ax[0], legend=False, zorder=4)
ax[0].set_ylim(-1.5, 6)
ax[0].set_yticks(np.arange(0, 6, 2))
ax[0].set_ylabel(r'Benthic $\delta^{18}$O')
ax[0].invert_yaxis()
ax[0].set_zorder(3)


#Omega (surface)
p2=[]
colours=['#FCB0B3', '#F93943', '#7EB2DD', '#445E93']
df_names=['CO2system_pH_omega65', 'CO2system_pH_alkalinity', 'CO2system_pH_dic', 'CO2system_pH_ccd']
df_labels=[r'$\Omega = 6.5$', 'Alkalinity-based', 'DIC-based', 'CCD-based']
for i, (p3_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p3_df].copy()
    temp_df['age']=temp_df['age']/1000
    p2.append(sns.lineplot(data=temp_df, x='age', y='saturation_state',ax=ax[1], 
                       color=c, label=lab, legend=False))
ax[1].set_ylim(3, 9)
ax[1].set_ylabel('Surface $\Omega_c$')



#Omega (benthic)
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p3=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='omega_logfit', 
                hue='ocean_basin', palette=custom_palette,ax=ax[2],zorder=3, s=20)
ax[2].set_ylabel(r'Deep $\Omega_c$')
#ax[1].set_ylim(25, 225)
ax[2].legend(fontsize=10, loc='upper left')
ax[2].set_ylim(0.2, 2)
ax[2].set_yticks([0.5, 1, 1.5, 2])



#CO3 (surface)
p4=[]
colours=['#FCB0B3', '#F93943', '#7EB2DD', '#445E93']
for i, (p5_df, lab, c) in enumerate(zip(df_names, df_labels, colours)):
    temp_df=Rae2022_dict[p5_df].copy()
    temp_df['age']=temp_df['age']/1000
    p4.append(sns.lineplot(data=temp_df, x='age', y='co3',ax=ax[3], 
                       color=c, label=lab))

colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p4b=sns.scatterplot(data=foram_df_B11_agecut, x='age_Ma', y='surface_CO3', 
                color='black', ax=ax[3], s=20, 
                legend=False, label='B/Ca-based')

ax[3].set_ylabel(r'Surface [CO$_3^{2-}$] ($\mu$mol/kg)')
ax[3].legend(loc='upper left', ncols=2)


#CO3 (benthic) converted
colours=['#1b9e77', '#d95f02', '#7570b3']
custom_palette=sns.set_palette(sns.color_palette(colours))
p5=sns.scatterplot(data=foram_df_B11, x='age_Ma', y='CO3', 
                hue='ocean_basin', palette=custom_palette, ax=ax[4], s=20, 
                legend=False, zorder=3)
#ax[3].legend(fontsize=10, loc='center left')
ax[4].set(ylabel='Deep [CO$_3^{2-}$] ($\mu$mol/kg)')
ax[4].set(xlabel='Age (Ma)')
ax[4].set_yticks([25, 50, 75, 100])
ax[4].set_ylim(5, 110)




ax[0].set_xlim(0, 70)
ax[4].set_xlabel('Age (Ma)')
ax[0].invert_xaxis()



#move the subplots together so that they slightly overlap
plt.subplots_adjust(hspace=-0.2)
#remove the box around the plots

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
    for epoch, bounds in epoch_boundaries.items():
        axis.axvline(bounds[0], color='lightgrey', linestyle='--', zorder=1)
#remove x ticks on all but the bottom plot

fig.savefig(fig_path/'EGU_fig11b.png', dpi=300)



## Z & T 2019

fig, ax = plt.subplots(figsize=(6, 6), nrows=1, sharex=True, constrained_layout=False)
sns.lineplot(data=ZandT['boron'], x='Age(Ma)', y='[CO32-](mumol/kg)', 
                ax=ax, color='black', label='Boron-based', markers=True, marker='o')
sns.lineplot(data=ZandT['CO3model'], x='age', y='CO3_CCD', ax=ax, color='red',  
                label='Geocarb III', markers=True, marker='o')
ax.set_xlim(0, 70)
ax.invert_xaxis()

ax2 = ax.twinx()  
ax2.plot(ZandT['boron']['Age(Ma)'].values, CO3_offset, color='blue', label='Offset', linestyle='--')
ax.set_xlabel('Age (Ma)')
ax.set_ylabel(r'[CO$_3^{2-}$] ($\mu$mol/kg)')
ax2.set_ylabel('Offset (mumol/kg)')

plt.title('Zeebe and Tyrrell (2019)')
ax2.legend(loc='lower right')
fig.savefig(fig_path/'ZandT2019.png', dpi=300)
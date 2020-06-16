'''
Script used for data analysis linked to Figure 2 and S2
'''

# %% General imports and setup

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plotnine import *
# IMPORTANT: use plotine v0.5.1, later versions might not work as expected
from scipy import stats
import scikit_posthocs
from collections import defaultdict
from scipy.spatial import distance
import seaborn as sns
from scipy.spatial.distance import cdist
from scipy import ndimage

# %% Figure 2D
# Data loading and plot

d_flat_mono_mt_ends = pd.read_csv('./data/microtubule_distribution_tips.csv')
p_microtubule_dist_ends = ggplot(d_flat_mono_mt_ends, aes('Length')) + geom_histogram(aes(y='..density..'),binwidth=0.7) + geom_density() + \
    theme_classic() + xlim(0,25)
print(p_microtubule_dist_ends)

# %% Figure 2E
# Data loading and plot

d_lgn_flat_rnai = pd.read_csv('./data/rnai_lgn_flat_monopolar.csv')

conditions=['control_rnai','lgn_rnai']
p_lgn_flat_rnai_treat_vel = ggplot(d_lgn_flat_rnai.loc[np.isin(d_lgn_flat_rnai['Treatment'], conditions)], aes('Treatment', 'Pole_velocity')) + \
    geom_boxplot(outlier_stroke=0, outlier_size=0) + geom_jitter(height=0) + theme_classic()
print(p_lgn_flat_rnai_treat_vel)

stats.mannwhitneyu(d_lgn_flat_rnai.loc[d_lgn_flat_rnai['Treatment']=='control_rnai','Pole_velocity'],
                  d_lgn_flat_rnai.loc[d_lgn_flat_rnai['Treatment']=='lgn_rnai','Pole_velocity'])


d_lgn_flat_noc = pd.read_csv('./data/flat_mono_nocodazole.csv')
conditions=['before_nocodazole','after_nocodazole']
p_lgn_flat_noc_vel = ggplot(d_lgn_flat_noc.loc[np.isin(d_lgn_flat_noc['Treatment'], conditions)], aes('Treatment', 'Centrosome_Velocity')) + \
    geom_boxplot(outlier_stroke=0, outlier_size=0) + geom_jitter(height=0) + theme_classic()
print(p_lgn_flat_noc_vel)

stats.mannwhitneyu(d_lgn_flat_noc.loc[d_lgn_flat_noc['Treatment']=='before_nocodazole','Centrosome_Velocity'],
                  d_lgn_flat_noc.loc[d_lgn_flat_noc['Treatment']=='after_nocodazole','Centrosome_Velocity'])

# %% Figure 2G
# Data loading and processing - This part of the code requires some computational time
d_flat_mono_lgn_imp_global = pd.read_csv('./data/data_lgn_importazole.csv')
#set(d_flat_mono_lgn_imp_global['cell_id'])

xy_res = 0.2369

d_flat_mono_lgn_imp = pd.DataFrame()
for cell_id in set(d_flat_mono_lgn_imp_global['cell_id']):
    print(str(cell_id)+' of '+str(max(d_flat_mono_lgn_imp_global['cell_id'])))
    d_flat_mono_lgn_imp_dna = pd.read_csv('./data/importazole/dna-id'+str(cell_id)+'.csv')
    d_flat_mono_lgn_imp_lgn = pd.read_csv('./data/importazole/lgn-id'+str(cell_id)+'.csv')

    d_flat_mono_lgn_imp_lgn['cell_id'] = cell_id
    d_flat_mono_lgn_imp_lgn['treatment'] = d_flat_mono_lgn_imp_global.loc[d_flat_mono_lgn_imp_global['cell_id']==cell_id,'treatment'].iloc[0]

    d_flat_mono_lgn_imp_lgn['min_dist'] = [np.min(cdist(d_flat_mono_lgn_imp_dna[['X','Y']],
                                        [x])) for x in np.asarray(d_flat_mono_lgn_imp_lgn[['X','Y']])]
    d_flat_mono_lgn_imp_lgn['min_dist'] *= xy_res

    d_flat_mono_lgn_imp_lgn['lgn'] = (d_flat_mono_lgn_imp_lgn['Value'] - \
                                      d_flat_mono_lgn_imp_global.loc[d_flat_mono_lgn_imp_global['cell_id']==cell_id,'lgn_background'].values) / \
                                     (d_flat_mono_lgn_imp_global.loc[d_flat_mono_lgn_imp_global['cell_id']==cell_id,'lgn_mean'].values - \
                                      d_flat_mono_lgn_imp_global.loc[d_flat_mono_lgn_imp_global['cell_id']==cell_id,'lgn_background'].values)

    d_flat_mono_lgn_imp = d_flat_mono_lgn_imp.append(d_flat_mono_lgn_imp_lgn)

d_flat_mono_lgn_imp['lgn_filt']=0
med_filt_rad = 5

for cell_id in set(d_flat_mono_lgn_imp['cell_id']):
    print(str(cell_id)+' of '+str(max(d_flat_mono_lgn_imp['cell_id'])))

    tmp_transformed_data = np.zeros((np.max(d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'X']),
                                     np.max(d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'Y'])))
    for i in range(0,len(d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id])):
        tmp_transformed_data[
            d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].X-1,
            d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].Y-1] = \
                d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].lgn

    tmp_filt_transformed_data = ndimage.median_filter(tmp_transformed_data, med_filt_rad)

    for i in range(0,len(d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id])):
        d_flat_mono_lgn_imp.loc[(d_flat_mono_lgn_imp.cell_id==cell_id)&\
                                (d_flat_mono_lgn_imp.X==d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].X)&\
                                (d_flat_mono_lgn_imp.Y==d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].Y), 'lgn_filt'] = \
            tmp_filt_transformed_data[d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].X-1,
                                  d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].iloc[i].Y-1]

d_flat_mono_lgn_imp['min_dist_round'] = np.round(d_flat_mono_lgn_imp['min_dist'],0)
g_d_flat_mono_lgn_imp = d_flat_mono_lgn_imp.groupby(['treatment','cell_id','min_dist_round']).mean().reset_index()

d_flat_mono_lgn_imp['is_high_lgn'] = False
tmp_quant = 0.975
for cell_id in set(d_flat_mono_lgn_imp.cell_id):
    #print(cell_id)
    tmp_thresh = d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].lgn_filt.quantile(tmp_quant)

    d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'is_high_lgn'] = \
        d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'lgn_filt'] >= tmp_thresh

g_d_flat_mono_lgn_imp_quant = d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.is_high_lgn].groupby(['treatment','cell_id']).min().reset_index()


# Plot for Figure 2G

p_dist_lgn_dist_hist_quant = ggplot(g_d_flat_mono_lgn_imp_quant.loc[g_d_flat_mono_lgn_imp_quant.is_high_lgn], aes('min_dist_round'))+ \
    geom_histogram(aes(y='..density..'),binwidth=2.5,center=1.25) +geom_density() + theme_classic() + facet_wrap('treatment') + \
    xlab('distance from DNA')
print(p_dist_lgn_dist_hist_quant)

# Statistical test for 2G
print(
    stats.mannwhitneyu(g_d_flat_mono_lgn_imp_quant.loc[g_d_flat_mono_lgn_imp_quant.treatment=='control','min_dist'],
                  g_d_flat_mono_lgn_imp_quant.loc[g_d_flat_mono_lgn_imp_quant.treatment=='importazole','min_dist'])
)





# Data preparation for Figure S2H

quant = np.arange(0.90,1.0,0.001)
lgn_range_mean = np.zeros(np.shape(quant))
quant_thresh = np.zeros(np.shape(quant))
i=0

for tmp_quant in quant:
    d_flat_mono_lgn_imp['is_high_lgn'] = False
    thresh = np.zeros(len(set(d_flat_mono_lgn_imp.cell_id)))
    j=0
    for cell_id in set(d_flat_mono_lgn_imp.cell_id):
        #print(cell_id)
        tmp_thresh = d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id].lgn_filt.quantile(tmp_quant)
        thresh[j]=tmp_thresh
        j+=1
        d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'is_high_lgn'] = \
            d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.cell_id==cell_id,'lgn_filt'] >= tmp_thresh

    g_d_flat_mono_lgn_imp_quant = d_flat_mono_lgn_imp.loc[d_flat_mono_lgn_imp.is_high_lgn].groupby(['treatment','cell_id']).min().reset_index()
    lgn_range_mean[i] = np.mean(g_d_flat_mono_lgn_imp_quant.loc[g_d_flat_mono_lgn_imp_quant.treatment=='control','min_dist'])
    quant_thresh[i] = np.mean(thresh)
    i+=1


d = {'quantile': quant, 'threshold':quant_thresh, 'lgn_inhibition_range': lgn_range_mean}
d_lgn_inh_range = pd.DataFrame(data=d)

# Plot of Figure S2H

p_quantile_lgn_inh_range = ggplot(d_lgn_inh_range, aes('quant','lgn_inhibition_range')) + geom_line() + theme_classic() +\
    ylim(0,8) + xlab('quantile') + ylab('lgn inhibition range')
print(p_quantile_lgn_inh_range)

# %% Figure S2G
# Data loading and plot

d_flat_mono_microtubules = pd.read_csv('./data/microtubule_distribution_intensity.csv')

for cell_id in set(d_flat_mono_microtubules.cell_id):
    max_tub_intensity = max(d_flat_mono_microtubules.loc[d_flat_mono_microtubules.cell_id==cell_id,'intensity'].values)
    d_flat_mono_microtubules.loc[d_flat_mono_microtubules.cell_id==cell_id,'tub_norm'] = \
        (d_flat_mono_microtubules.loc[d_flat_mono_microtubules.cell_id==cell_id,'intensity'] - d_flat_mono_microtubules.loc[d_flat_mono_microtubules.cell_id==cell_id,'background']) / \
        (max_tub_intensity - d_flat_mono_microtubules.loc[d_flat_mono_microtubules.cell_id==cell_id,'background'])

p_microtubule_dist_intensity = ggplot(
    d_flat_mono_microtubules,
    aes('distance', 'tub_norm')
) + geom_line(aes(colour='factor(cell_id)')) + stat_summary() + theme_classic() + xlim(0,15)
print(p_microtubule_dist_intensity)

'''
Script used for data analysis linked to Figure 1 and S1
'''

# %% General imports and setup

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plotnine import *
# IMPORTANT: use plotine v0.5.1, later versions might not work as expected
from scipy import stats
#import scikit_posthocs
from collections import defaultdict
from scipy.spatial import distance
import seaborn as sns

initial_cell_shape_time = -9 # minutes
spindle_formation_time = 9 # minutes
metaphase_time = 21 # minutes

# %% Figure 1C
# Data loading and pre-processing

d_flat_mono = pd.read_csv('./data/flat_monopolar.csv')

d_flat_mono['initial_cell_x'] = np.nan
d_flat_mono['initial_cell_y'] = np.nan
d_flat_mono['initial_cell_angle'] = np.nan

for cell_id in set(d_flat_mono['cell_id']):
    d_flat_mono.loc[
        (d_flat_mono['cell_id']==cell_id),
        'initial_cell_x'
    ] = d_flat_mono.loc[
            (d_flat_mono['time']==initial_cell_shape_time) &
            (d_flat_mono['cell_id']==cell_id),'cell_x'].values

for cell_id in set(d_flat_mono['cell_id']):
    d_flat_mono.loc[
        (d_flat_mono['cell_id']==cell_id),
        'initial_cell_y'
    ] = d_flat_mono.loc[
            (d_flat_mono['time']==initial_cell_shape_time) &
            (d_flat_mono['cell_id']==cell_id),'cell_y'].values

for cell_id in set(d_flat_mono['cell_id']):
    d_flat_mono.loc[
        (d_flat_mono['cell_id']==cell_id),
        'initial_cell_angle'
    ] = d_flat_mono.loc[
            (d_flat_mono['time']==initial_cell_shape_time) &
            (d_flat_mono['cell_id']==cell_id),'cell_angle'].values

d_flat_mono['centrosomes_x'] = d_flat_mono['centrosome_mid_x'] - d_flat_mono['initial_cell_x']
d_flat_mono['centrosomes_y'] = -(d_flat_mono['centrosome_mid_y'] - d_flat_mono['initial_cell_y'])
d_flat_mono['centrosomes_x_norm'] = d_flat_mono['centrosomes_x'] * np.cos(d_flat_mono['initial_cell_angle']/180*np.pi) + \
                                    d_flat_mono['centrosomes_y'] * np.sin(d_flat_mono['initial_cell_angle']/180*np.pi)
d_flat_mono['centrosomes_y_norm'] = - d_flat_mono['centrosomes_x'] * np.sin(d_flat_mono['initial_cell_angle']/180*np.pi) + \
                                    d_flat_mono['centrosomes_y'] * np.cos(d_flat_mono['initial_cell_angle']/180*np.pi)
d_flat_mono['centrosome_dna_dist'] = ((d_flat_mono.centrosome_mid_x - d_flat_mono.dna_x)**2 + (d_flat_mono.centrosome_mid_y - d_flat_mono.dna_y)**2)**0.5

# Plot of figure 1C

p_centrosomes_trajectory = ggplot(d_flat_mono.loc[(d_flat_mono['time']>0)&(d_flat_mono['time']<=72)], \
    aes('centrosomes_x_norm', 'centrosomes_y_norm', group='cell_id', colour='factor(cell_id)')) + \
    geom_path() + coord_fixed() + xlim(-60,60) + ylim(-40,40) + theme_classic()
print(p_centrosomes_trajectory)

# %% Figure 1G
# Data loading

d_fn_peg = pd.read_csv('./data/pdms_confinement_fn_peg.csv')

pdms_fn = d_fn_peg.loc[d_fn_peg['Treatment']=='pdms_fn',['Cell','Treatment','Velocity']]
pdms_peg = d_fn_peg.loc[d_fn_peg['Treatment']=='pdms_peg',['Cell','Treatment','Velocity']]
rap1_fn = d_flat_mono.loc[d_flat_mono['time']>spindle_formation_time, ['cell_id','condition','centrosome_mid_vel']]
rap1_fn.cell_id += 10
rap1_fn.columns = pdms_fn.columns
d_rap1_pdms_fn_peg = pdms_fn.append(pdms_peg, ignore_index=True).append(rap1_fn, ignore_index=True)

# Plot of figure 1G
g_d_rap1_pdms_fn_peg = d_rap1_pdms_fn_peg.groupby(['Treatment','Cell']).mean().reset_index()
p_grouped_fn_peg_velocity = ggplot(g_d_rap1_pdms_fn_peg, aes('Treatment', 'Velocity')) + geom_boxplot(outlier_size=0,
    outlier_stroke=0) + geom_jitter(height=0) + ylim(0,2.5)+ theme_classic()
print(p_grouped_fn_peg_velocity)

# Stats for figure 1G

s_kru_g = stats.kruskal(g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='pdms_fn','Velocity'],
              g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='pdms_peg','Velocity'],
              g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='rap1_stlc','Velocity'], nan_policy='omit')

s_post_g = scikit_posthocs.posthoc_conover([g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='pdms_fn','Velocity'],
                                 g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='pdms_peg','Velocity'],
                                 g_d_rap1_pdms_fn_peg.loc[g_d_rap1_pdms_fn_peg.Treatment=='rap1_stlc','Velocity']])

print(s_kru_g, s_post_g)

# %% Figure 1I
# Data loading
d_flat_mono_lines = pd.read_csv('./data/flat_monopolar_spindle_line_patterns.csv')

# Plot of figure 1I
p_cells_on_lines_position_velocity = ggplot(d_flat_mono_lines, \
    aes('x_pole', 'v_pole', \
    group='id', colour='factor(id)')) +\
    geom_path() + theme_classic()
print(p_cells_on_lines_position_velocity)

# %% Figure S1B
# Data loading
data_round_bipolar = pd.read_csv('./data/rounding_cells_bipolar_spindles.csv')
spindle_final_shape = 15  ## minutes

# Plot of figure S1B
p_round_bipolar_dna_major = ggplot(data_round_bipolar.loc[(data_round_bipolar['time']>spindle_final_shape)], aes(
    'dna_length')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) +xlim(0,35) + theme_classic()
p_round_bipolar_dna_minor = ggplot(data_round_bipolar.loc[(data_round_bipolar['time']>spindle_final_shape)], aes(
    'dna_width')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) +xlim(0,35) +theme_classic()
p_round_bipolar_centrosome_dist = ggplot(data_round_bipolar.loc[(data_round_bipolar['time']>spindle_final_shape)], aes(
    'spindle_length')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) +xlim(0,20) +theme_classic()

print(p_round_bipolar_dna_major,
      p_round_bipolar_dna_minor,
      p_round_bipolar_centrosome_dist)

# %% Figure S1D
# Plot of figure S1D

p_monopolar_geometry_dna_major = ggplot(d_flat_mono.loc[(d_flat_mono['time']>spindle_final_shape)], aes(
    'dna_length')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) + xlim(0,35)+ theme_classic()
p_monopolar_geometry_dna_minor = ggplot(d_flat_mono.loc[(d_flat_mono['time']>spindle_final_shape)], aes(
    'dna_width')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) + xlim(0,35) +theme_classic()
p_monopolar_geometry_centrosome_dist = ggplot(d_flat_mono.loc[(d_flat_mono['time']>spindle_final_shape)], aes(
    'centrosome_distance')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) + xlim(0,20) +theme_classic()
p_monopolar_geometry_pole_dna_dist = ggplot(d_flat_mono.loc[(d_flat_mono['time']>spindle_final_shape)],
    aes('centrosome_dna_dist')) + geom_histogram(aes(y='..density..'),binwidth=2,center=1) + xlim(0,20) +theme_classic()

print(p_monopolar_geometry_dna_major,p_monopolar_geometry_dna_minor,
      p_monopolar_geometry_centrosome_dist, p_monopolar_geometry_pole_dna_dist)

# %% Figure S1F-G-H-I
# Data loading for S1F-G-H
dVel = pd.read_csv("./data/pole-vel-4-conditions-grouped.csv")
cols = ['condition', 'id', 'time', 'spindle_pole_a_vel', 'spindle_pole_b_vel']
dVel = dVel[cols]
dVel['pole-vel'] = np.mean(dVel[['spindle_pole_a_vel','spindle_pole_b_vel']], axis=1)

# Plots of figure S1F
pVelControl = sns.distplot(dVel.loc[dVel['condition']=='control', "pole-vel"])
pVelControl.set(ylim=(0,2),xlim=(0,5))

pVelRap1Stlc = sns.distplot(dVel.loc[dVel['condition']=='rap1_stlc', "pole-vel"])
pVelRap1Stlc.set(ylim=(0,2),xlim=(0,5))

pVelRap1 = sns.distplot(dVel.loc[dVel['condition']=='rap1', "pole-vel"])
pVelRap1.set(ylim=(0,2),xlim=(0,5))

pVelStlc = sns.distplot(dVel.loc[dVel['condition']=='stlc', "pole-vel"])
pVelStlc.set(ylim=(0,2),xlim=(0,5))

# Plots of figure S1G

groupedCondIdMean = dVel[dVel['time']>spindle_final_shape].groupby(['condition','id']).mean().reset_index()
pPoleVelocity = sns.catplot(x='condition', y='pole-vel', kind='box',
                            order=['control', 'rap1', 'stlc', 'rap1_stlc'],
                            color='w',
                            data=groupedCondIdMean).set(ylim=(0,2.5))
sns.swarmplot(x='condition', y='pole-vel',
        order=['control', 'rap1', 'stlc', 'rap1_stlc'],
        color='k',
        data=groupedCondIdMean, ax=pPoleVelocity.ax)

# Plots of figure S1H

groupedCondIdMax = dVel[dVel['time']>spindle_final_shape].groupby(['condition','id']).max().reset_index()
pPoleVelocityMax = sns.catplot(x='condition', y='pole-vel', kind='box',
                            order=['control', 'rap1', 'stlc', 'rap1_stlc'],
                            color='w',
                            data=groupedCondIdMax).set(ylim=(0,5.2))
sns.swarmplot(x='condition', y='pole-vel',
        order=['control', 'rap1', 'stlc', 'rap1_stlc'],
        color='k',
        data=groupedCondIdMax, ax=pPoleVelocityMax.ax)

# Data loading for S1I

dPos = pd.read_csv("./data/pole-positions-4-conditions.csv")
dPos.columns = ['raw','cell','condition','id','time','mid-x','mid-y']

dMaxDistDict = defaultdict(list)
for i in set(dPos['id']) - set([38]):
    meanX = dPos.loc[dPos['id'] == i, 'mid-x'].mean()
    meanY = dPos.loc[dPos['id'] == i, 'mid-y'].mean()
    maxDist = np.array([distance.euclidean(el, (meanX, meanY))
            for el in list(zip(
            dPos.loc[(dPos['id'] == i) & (dPos['time'] >= 15), 'mid-x'],
            dPos.loc[(dPos['id'] == i) & (dPos['time'] >= 15), 'mid-y']
            ))]).max()
    dMaxDistDict['id'].append(i)
    dMaxDistDict['condition'].append(dPos.loc[dPos['id'] == i, 'condition'].iloc[0])
    dMaxDistDict['max-dist'].append(maxDist)

# Plot of figure S1I
dMaxDist = pd.DataFrame(dMaxDistDict)
pMaxDist = sns.catplot(x='condition', y='max-dist', kind='box',
                            order=['control', 'stlc', 'rap1', 'rap1_stlc'],
                            color='w',
                            data=dMaxDist)
sns.swarmplot(x='condition', y='max-dist',
        order=['control', 'stlc', 'rap1', 'rap1_stlc'],
        color='k',
        data=dMaxDist, ax=pMaxDist.ax)

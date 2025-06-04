#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 19:34:40 2025

@author: brendonmcguinness
"""
import cometspy as c
from cobra.io import load_model, read_sbml_model
import numpy as np
import matplotlib.pyplot as plt
from helpers import nicheOverlapSaveTimeSeriesMultipleChooseTime, returnSCMatrix

exchange_rxns = ["EX_glc__D_e", "EX_xyl__D_e"]
exchange_rxns = ["EX_glc__D_e", "EX_succ_e"]
#exchange_rxns = ["EX_glc__D_e", "EX_fru_e"]



# concentrations: both at 0.1
s1_conc = np.array([0.1])
s2_conc = np.array([0.1])
mut1 = c.model(read_sbml_model('ecoli_k12_mg1655.xml'))
#bs = c.model(read_sbml_model('bacillus_subtilis.xml'))
mut2 = c.model(load_model("iYO844"))
# cross-bound non-target rxns to zero
ko1=-5
ko2=-10
for rxn in exchange_rxns:
    mut1.change_bounds(rxn, ko1, 1000)
    mut2.change_bounds(rxn, ko2, 1000)
 

source1 = exchange_rxns[0]
source2 = exchange_rxns[1] 
source_list = [source1,source2]
#mut1.change_bounds(exchange_rxns[0], ko1, 1000)
#mut2.change_bounds(exchange_rxns[1], ko2, 1000)    
                                                                  
mut1.id = "ecoli"
mut2.id = "bsubtilis"
# apply bounds
mut1.change_bounds(source1, ko1, 1000)
mut2.change_bounds(source2, ko2, 1000)

M = 11
s_wt, s_mut, wt_freq, media_list = returnSCMatrix(mut1,mut2,source1,source2,s1_conc,s2_conc,N=M,init_ratio_diff = 1e-3, init_mean_pop=5e-6)

R1 = []
R2 = []
for media in media_list:
    last = media[media.cycle == media.cycle.max()]
    met1 = source1[3:]
    met2 = source2[3:]
    
    # strain-1’s resource
    mask1 = last['metabolite'] == met1
    if mask1.any():
        R1_val = last.loc[mask1, 'conc_mmol'].iloc[0]
    else:
        R1_val = 0.0   # completely consumed
    
    # strain-2’s resource
    mask2 = last['metabolite'] == met2
    if mask2.any():
        R2_val = last.loc[mask2, 'conc_mmol'].iloc[0]
    else:
        R2_val = 0.0
    
    R1.append(R1_val)
    R2.append(R2_val)

A_i = np.zeros((M,len(source_list)))
for idx in range(len(wt_freq)):
        for k,s in enumerate(source_list):
            #A_i[idx] = np.trapz(media_list[idx].loc[media_list[idx].metabolite==source1[3:], 'conc_mmol'],
            #              x=media_list[idx].cycle * 0.1)
            sub = media_list[idx][media_list[idx].metabolite == s[3:]]
            t = sub['cycle'].to_numpy() * 0.1   # now t.shape == sub.shape[0]
            c = sub['conc_mmol'].to_numpy()
            A_i[idx,k] = np.trapz(c, x=t)

NO = (s_wt[0]-s_wt[-1])/(wt_freq[0]-wt_freq[-1])
ND = abs(NO)
coex_strength = np.min(np.array([s_wt[0],s_mut[1]]))
FD = s_wt[0]-s_mut[1]
sc_wt_flat = np.ones(M) * s_wt[0] #np.array([s_wt[0],s_wt[0]])
sc_mut_flat = np.ones(M) * s_mut[-1] #np.array([s_mut[1],s_mut[1]])

plt.figure()
plt.scatter(wt_freq,s_wt,label='strain 1')
plt.plot(wt_freq,s_wt,color='tab:blue')

plt.scatter(wt_freq[::-1],s_mut,label='strain 2')
plt.plot(wt_freq[::-1],s_mut,color='tab:orange')

plt.plot(wt_freq,sc_wt_flat,color='tab:blue',linestyle='dashed',label='SC when rare (strain 1)')
plt.plot(wt_freq[::-1],sc_mut_flat,color='tab:orange',linestyle='dashed',label='SC when rare (strain 2)')

#plt.axhline(s_wt[0],color='k',linestyle='dashed')
#plt.axhline(s_mut[1],color='k',linestyle='dashed')

plt.annotate('Niche difference='+f"{ND:.3f}", xy=(349, 250), xycoords='axes points',
            size=14, ha='right', va='top',
            bbox=dict(boxstyle='round', fc='w'))
#plt.semilogx([0.01,0.5,0.99],[0.006399633,0.01216757,0.032638824])
#plt.legend(['Niche overlap='+f"{niche_overlap:.3f}"])
#plt.ylim(-1,1)
#plt.xscale('log')
#plt.yscale('log')
plt.ylabel('selection coefficient')
plt.xlabel('Initial relative abundance (Strain 1:Strain 2)')
plt.legend()
plt.show()

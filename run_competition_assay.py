#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: brendonmcguinness
"""
#monoculture

# Start by loading required packages, including the COMETS toolbox
import cometspy as c
#from cobra import test
import cobra
from cobra.io import load_model,read_sbml_model
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from cobra.medium import minimal_medium
from helpers import flux_cosim

#os.environ['GUROBI_HOME'] = 'C:\\gurobi902\\win64'
os.environ['COMETS_HOME'] = '/Applications/COMETS'
os.environ['GUROBI_COMETS_HOME'] = '/Library/gurobi1003/macos_universal2'
os.environ['GRB_LICENSE_FILE']='/Library/gurobi1003/macos_universal2/gurobi.lic' # load the models and perform the mutation
#
#wt = c.model(test.create_test_model("ecoli"))
#ecoli_model = load_model("iJO1366")
ec = c.model(read_sbml_model('ecoli_k12_mg1655.xml'))
#bs = c.model(read_sbml_model('bacillus_subtilis.xml'))
ec.id = 'ecoli'
X = -10 


ec.change_bounds('EX_succ_e', -20,1000)
ec.change_bounds('EX_glc__D_e',X,1000)

bs = c.model(load_model("iYO844"))
bs.id = 'bsubtilis'
bs.change_bounds('EX_succ_e', -X,1000)
bs.change_bounds('EX_glc__D_e',-20,1000)


# set its initial biomass, 5e-6 gr at coordinate [0,0]
ec.initial_pop = [0, 0, 5e-6]
bs.initial_pop = [0, 0, 5e-6]

# create an empty layout
test_tube = c.layout()

# add the models to the test tube
test_tube.add_model(ec)
test_tube.add_model(bs)


# Add glucose to the media 
glc_molarmass = 180.156
succ_molarmass = 118.09
glc_conc = 0.1 # mol/L
succ_conc = 0.1 # mol/L
glc_gram = glc_conc*glc_molarmass*0.001
succ_gram = succ_conc*succ_molarmass*0.001 #1mL=1cm^3

test_tube.set_specific_metabolite('succ_e', succ_conc)
test_tube.set_specific_metabolite('glc__D_e', glc_conc)
#test_tube.set_specific_metabolite('lac__D_e', glc_conc)
#test_tube.set_specific_metabolite('', glc_conc)



# Add typical trace metabolites and oxygen coli as static

trace_metabolites = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e', 'h_e', 'k_e', 'h2o_e', 'mg2_e',
                     'mn2_e', 'mobd_e', 'na1_e', 'ni2_e', 'nh4_e', 'o2_e', 'pi_e', 'so4_e', 'zn2_e']
max_time = 240
for i in trace_metabolites:
    test_tube.set_specific_metabolite(i, 1000)
    test_tube.set_specific_static(i, 1000)
    
comp_params = c.params()
comp_params.set_param('maxCycles', max_time)
comp_params.set_param('writeMediaLog', True)
#comp_params.set_param('writeFluxLog', True)
comp_params.set_param('FluxLogRate',1)
comp_params.set_param('MediaLogRate',1)

comp_assay = c.comets(test_tube, comp_params)
comp_assay.run()

biomass = comp_assay.total_biomass # * 1000 / 0.05
biomass['t'] = biomass['cycle'] * comp_assay.parameters.all_params['timeStep']


myplot = biomass.drop(columns=['cycle']).plot(x = 't')
myplot.set_ylabel("Biomass (gr.)")
myplot.set_yscale('log')
myplot.legend([ec.id,bs.id])

m_wt = np.log(biomass.ecoli[max_time]/biomass.ecoli[0])
y_glc = biomass.ecoli[max_time] / glc_gram
print('Yield=',y_glc)
#m_mut = np.log(biomass.succ_KO[240]/biomass.succ_KO[0])

media = comp_assay.media.copy()
media = media[media.conc_mmol<900]

glc_ts = media[media['metabolite']=='glc__D_e']
succ_ts = media[media['metabolite']=='succ_e']
ac_ts = media[media['metabolite']=='ac_e']
carbons = pd.concat([glc_ts, succ_ts, ac_ts])

fig, ax = plt.subplots()
media.groupby('metabolite',sort=False).plot(x='cycle', ax =ax, y='conc_mmol')
ax.set_xlim((0,max_time))
#ax.legend(('glucose','succinate','acetate'))
ax.legend(media.metabolite.unique())
ax.set_ylabel("Concentration (mmol)")
"""
wt_df = comp_assay.fluxes_by_species['wt']
ex_glc = wt_df.EX_glc__D_e
ex_succ = wt_df.EX_succ_e
plt.figure()
plt.plot(wt_df.cycle/10,np.cumsum(wt_df.EX_co2_e),label='co2')
plt.plot(wt_df.cycle/10,np.cumsum(wt_df.EX_ac_e),label='acetate')
plt.ylabel('cumulative flux')
plt.xlabel('time (hours)')
plt.legend()
plt.show()
"""
#plt.plot(flux_cosim(wt,mut,comp_assay))


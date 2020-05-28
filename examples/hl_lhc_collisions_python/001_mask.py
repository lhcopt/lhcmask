import os

#from cpymad.madx import Madx
from madxp import Madxp as Madx
import pymasktools as pmt

beam_to_configure = 1

pmt.make_links(force=True, links_dict={
    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
    'modules': 'tracking_tools/modules',
    'tools': 'tracking_tools/tools',
    'beambeam_macros': 'tracking_tools/beambeam_macros',
    'errors': 'tracking_tools/errors'})


mad = Madx()

# Build sequence (to become python function)
mad.input(f'mylhcbeam = {beam_to_configure}')
mad.call('hl14_thin.madx')

# Make optics (to become python function)
mad.call('hl14_collision_optics.madx')

# Load parameters 
from parameters import parameters
pmt.checks_on_parameter_dict(parameters)


mad.set_variables_from_dict(params=parameters)

# Call mask modules
mad.call("modules/submodule_01a_preparation.madx")
mad.call("modules/submodule_01b_beam.madx")

var_dict = mad.get_variables_dicts()

import optics_specific_tools as ost
sequences_to_check = ['lhcb1', 'lhcb2']
twiss_fname = 'twiss_from_optics'
twiss_dfs = {}
for ss in sequences_to_check:
    mad.use(ss)
    mad.twiss()
    tdf = mad.get_table_df('twiss')
    twiss_dfs[ss] = tdf

for ss in sequences_to_check:
    tt = twiss_dfs[ss]
    if twiss_fname is not None:
        tt.to_parquet(twiss_fname + f'_{ss}.parquet')


for ss in sequences_to_check:
    tt = twiss_dfs[ss]
    ost.check_beta_at_ips_against_madvars(beam=ss[-1],
            twiss_df=tt,
            variable_dicts=var_dict,
            tol=1e-3)
print('IP beta test against knobs passed!')

separations_to_check = [
        {'element_name': 'ip1:1', 'plane':'v', 'knob': 'on_sep1', 'tol': 1e-3}]


prrr

mad.call("modules/submodule_01c_phase.madx")



#mad.call("modules/module_02_lumilevel.madx")
#mad.call("modules/module_03_beambeam.madx")
#mad.call("modules/module_04_errors.madx")
#mad.call("modules/module_05_tuning.madx")
#mad.call("modules/module_06_generate.madx")

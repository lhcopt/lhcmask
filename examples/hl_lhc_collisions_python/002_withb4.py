import os

import numpy as np

#from cpymad.madx import Madx
from madxp import Madxp as Madx
import pymasktools as pmt
import optics_specific_tools as ost

beam_to_configure = 1
optics_file = 'hl14_collision_optics.madx' #15 cm
save_intermediate_twiss = True
sequences_to_check = ['lhcb1', 'lhcb2']

# Tolarances for checks [ip1, ip2, ip5, ip8]
tol_beta = [1e-3, 10e-2, 1e-3, 1e-2]
tol_sep = [1e-6, 1e-6, 1e-6, 1e-6]

pmt.make_links(force=True, links_dict={
    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
    'modules': 'tracking_tools/modules',
    'tools': 'tracking_tools/tools',
    'beambeam_macros': 'tracking_tools/beambeam_macros',
    'errors': 'tracking_tools/errors'})


mad = Madx()
# Build sequenc
ost.build_sequence(mad, beam=beam_to_configure)
# Apply optics
ost.apply_optics(mad, optics_file=optics_file)

# Check and load parameters 
from parameters import parameters
pmt.checks_on_parameter_dict(parameters)
mad.set_variables_from_dict(params=parameters)

# Prepare sequences and attach beam
mad.call("modules/submodule_01a_preparation.madx")
mad.call("modules/submodule_01b_beam.madx")

# Test machine before any change
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_from_optics',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=True, check_separations_at_ips=True)

# Set phase, apply and save crossing
mad.call("modules/submodule_01c_phase.madx")
mad.call("modules/submodule_01d_crossing.madx")

# Test flat machine
mad.input('exec, crossing_disable')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_no_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=True, check_separations_at_ips=True)
# Check flatness
flat_tol = 1e-6
for ss in twiss_dfs.keys():
    tt = twiss_dfs[ss]
    assert np.max(np.abs(tt.x)) < flat_tol
    assert np.max(np.abs(tt.y)) < flat_tol

# Check machine after crossing restore
mad.input('exec, crossing_restore')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_with_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=True, check_separations_at_ips=True)

# Call other modules
mad.use(f'lhcb{beam_to_configure}')
mad.call("modules/module_02_lumilevel.madx")

# Generate b4
mad_b2 = mad
mad_b4 = Madx()
ost.build_sequence(mad_b4, beam=4)

var_dicts_b2 = mad_b2.get_variables_dicts()
var_dicts_b4 = mad_b4.get_variables_dicts()


b2_const=var_dicts_b2['constants']
b4_const=var_dicts_b4['constants']
for nn in b2_const.keys():
    if nn[0]=='_':
        print(f'The constant {nn} cannot be assigned!')
    else:
        if nn not in b4_const.keys():
            mad_b4.input(f'const {nn}={b2_const[nn]}')

# %% INDEPENDENT
b2_indep=var_dicts_b2['independent_variables']
b4_indep=var_dicts_b4['independent_variables']
for nn in b2_indep.keys():
    mad_b4.input(f'{nn}={b2_indep[nn]}')

# %% DEPENDENT
b2_dep=var_dicts_b2['dependent_variables_expr']
b4_dep=var_dicts_b4['dependent_variables_expr']
for nn in b2_dep.keys():
    mad_b4.input(f'{nn}:={str(b2_dep[nn])}')

# bv_aux and my my lhcbeam need to be defined explicitly
mad_b4.input(f'bv_aux=-1')
mad_b4.input(f'mylhcbeam=4')

# Attach beam
mad_b4.input(str(mad.sequence['lhcb2'].beam))
mad_b4.use('lhcb2')
mad_b4.sequence['lhcb2'].beam['bv']=1

# %% CHECKS
var_dicts_b2 = mad_b2.get_variables_dicts()
var_dicts_b4 = mad_b4.get_variables_dicts()

b2_const=var_dicts_b2['constants']
b4_const=var_dicts_b4['constants']
for nn in b4_const.keys():
    assert b2_const[nn] == b4_const[nn]

for nn in b2_const.keys():
    if nn not in b4_const.keys():
        print(f'Warning: b2 const {nn}={b2_const[nn]} is not in b4.')

b2_indep=var_dicts_b2['independent_variables']
b4_indep=var_dicts_b4['independent_variables']
for nn in b2_indep.keys():
    if str(nn) in 'bv_aux mylhcbeam':
        continue
    assert b4_indep[nn] == b2_indep[nn]

for nn in b4_indep.keys():
    if nn not in b2_indep.keys():
        print(f'Warning: b4 indep {nn}={b4_indep[nn]} is not in b2.')

b2_dep=var_dicts_b2['dependent_variables_expr']
b4_dep=var_dicts_b4['dependent_variables_expr']
for nn in b2_dep.keys():
    if str(nn) in 'bv_aux mylhcbeam':
        continue
    assert str(b4_dep[nn]) == str(b2_dep[nn])

for nn in b4_dep.keys():
    if nn not in b2_dep.keys():
        print(f'Warning: b4 dep {nn}={str(b4_dep[nn])} is not in b2.')



prrrrr

mad.call("modules/module_03_beambeam.madx")
mad.call("modules/module_04_errors.madx")
mad.call("modules/module_05_tuning.madx")
mad.call("modules/module_06_generate.madx")

import os

import numpy as np

#from cpymad.madx import Madx
from madxp import Madxp as Madx
import pymasktools as pmt
import optics_specific_tools as ost

# Select mode
#mode = 'b1_without_bb'
mode = 'b1_with_bb'
#mode = 'b1_with_bb_legacy_macros'
#mode = 'b4_without_bb'
#mode = 'b4_without_bb_from_b2'
#mode = 'b4_with_bb_from_b2'


optics_file = 'hl14_collision_optics.madx' #15 cm

check_betas_at_ips = True
check_separations_at_ips = True
save_intermediate_twiss = True

if mode=='b1_without_bb':
    beam_to_configure = 1
    sequences_to_check = ['lhcb1', 'lhcb2']
    generate_b4_from_b2 = False
    include_bb = False
    bb_legacy_macros = False
elif mode=='b1_with_bb':
    beam_to_configure = 1
    sequences_to_check = ['lhcb1', 'lhcb2']
    generate_b4_from_b2 = False
    include_bb = True
    bb_legacy_macros = False
elif mode=='b1_with_bb_legacy_macros':
    beam_to_configure = 1
    sequences_to_check = ['lhcb1', 'lhcb2']
    generate_b4_from_b2 = False
    include_bb = True
    bb_legacy_macros = True
elif mode == 'b4_without_bb':
    beam_to_configure = 4
    sequences_to_check = ['lhcb2']
    generate_b4_from_b2 = False
    include_bb = False
    bb_legacy_macros = False
    check_separations_at_ips = False
elif mode == 'b4_without_bb_from_b2':
    beam_to_configure = 1
    sequences_to_check = ['lhcb1', 'lhcb2']
    generate_b4_from_b2 = True
    include_bb = False
    bb_legacy_macros = False
elif mode == 'b4_with_bb_from_b2':
    beam_to_configure = 1
    sequences_to_check = ['lhcb1', 'lhcb2']
    generate_b4_from_b2 = True
    include_bb = False
    bb_legacy_macros = False





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

# Build sequence
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
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)

# Set phase, apply and save crossing
mad.call("modules/submodule_01c_phase.madx")
mad.call("modules/submodule_01d_crossing.madx")

# Test flat machine
mad.input('exec, crossing_disable')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_no_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)
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
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)

# Call leveling module
mad.use(f'lhcb{beam_to_configure}')
mad.call("modules/module_02_lumilevel.madx")

# Generate b4
if generate_b4_from_b2:
    mad_b2 = mad
    mad_b4 = Madx()
    ost.build_sequence(mad_b4, beam=4)
    pmt.configure_b4_from_b2(mad_b4, mad_b2)

    twiss_dfs_b2, other_data_b2 = ost.twiss_and_check(mad_b2,
            sequences_to_check=['lhcb2'],
            tol_beta=tol_beta, tol_sep=tol_sep,
            twiss_fname='twiss_b2_for_b4check',
            save_twiss_files=save_intermediate_twiss,
            check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

    twiss_dfs_b4, other_data_b4 = ost.twiss_and_check(mad_b4,
            sequences_to_check=['lhcb2'],
            tol_beta=tol_beta, tol_sep=tol_sep,
            twiss_fname='twiss_b4_for_b4check',
            save_twiss_files=save_intermediate_twiss,
            check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)


# Install beam beam
import beambeam as bb
bb_dfs = bb.generate_bb_dataframes(mad,
    ip_names=['ip1', 'ip2', 'ip5', 'ip8'],
    numberOfLRPerIRSide=[25, 20, 25, 20],
    harmonic_number=35640,
    bunch_spacing_buckets=10,
    numberOfHOSlices=11,
    bunch_population_ppb=None,
    sigmaz_m=None)

prrrrr




# mad.call("modules/module_03_beambeam.madx")
mad.call("modules/module_04_errors.madx")
mad.call("modules/module_05_tuning.madx")
mad.call("modules/module_06_generate.madx")

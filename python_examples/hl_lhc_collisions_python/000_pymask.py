import os
import sys
import pickle

import numpy as np

from config import python_parameters, mask_parameters, knob_parameters


#####################################################
# Read general configurations and setup envirnoment #
#####################################################

mode = python_parameters['mode']
tol_beta = python_parameters['tol_beta']
tol_sep = python_parameters['tol_sep']
flat_tol = python_parameters['tol_co_flatness']
links = python_parameters['links']
optics_file = python_parameters['optics_file']
check_betas_at_ips = python_parameters['check_betas_at_ips']
check_separations_at_ips = python_parameters['check_separations_at_ips']
save_intermediate_twiss = python_parameters['save_intermediate_twiss']
force_leveling= python_parameters['force_leveling']

# Make links
for kk in links.keys():
    if os.path.exists(kk):
        os.remove(kk)
    os.symlink(os.path.abspath(links[kk]), kk)

# Execute customization script if present
os.system('bash customization.bash')

# Import pymask
sys.path.append('./modules')
import pymask as pm
import pymask.luminosity as lumi

# Import user-defined optics-specific tools
import optics_specific_tools as ost

######################################
# Check parameters and activate mode #
######################################

# Check and load parameters 
pm.checks_on_parameter_dict(mask_parameters)

# Define configuration
(beam_to_configure, sequences_to_check, sequence_to_track, generate_b4_from_b2,
    track_from_b4_mad_instance, enable_bb_python, enable_bb_legacy,
    force_disable_check_separations_at_ips,
    ) = pm.get_pymask_configuration(mode)

if force_disable_check_separations_at_ips:
    check_separations_at_ips = False

if not(enable_bb_legacy) and not(enable_bb_python):
    mask_parameters['par_on_bb_switch'] = 0.


########################
# Build MAD-X instance #
########################

# Start mad
Madx = pm.Madxp
mad = Madx()

# Build sequence
ost.build_sequence(mad, beam=beam_to_configure)

# Apply optics
ost.apply_optics(mad, optics_file=optics_file)

# Pass parameters to mad
mad.set_variables_from_dict(params=mask_parameters)

# Prepare auxiliary mad variables
mad.call("modules/submodule_01a_preparation.madx")

# Attach beams to sequences
mad.call("modules/submodule_01b_beam.madx")

# Test machine before any change
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_from_optics',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips,
        check_separations_at_ips=check_separations_at_ips)

# Set IP1-IP5 phase and store corresponding reference
mad.call("modules/submodule_01c_phase.madx")

# Set optics-specific knobs
ost.set_optics_specific_knobs(mad, knob_parameters, mode)

# Crossing-save and some reference measurements
mad.input('exec, crossing_save')
mad.call("modules/submodule_01e_final.madx")


#################################
# Check bahavior of orbit knobs #
#################################

# Check flat machine
mad.input('exec, crossing_disable')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_no_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)

# Check orbit flatness
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


#################################
# Set luminosity in IP2 and IP8 #
#################################


if enable_bb_legacy or mode=='b4_without_bb':
    mad.use(f'lhcb{beam_to_configure}')
    if mode=='b4_without_bb':
        print('Leveling not working in this mode!')
    else:
        # Luminosity levelling
        print('Luminosities before leveling (crab cavities are not considered):')
        lumi.print_luminosity(mad, twiss_dfs,
                mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
                mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])

        mad.call("modules/module_02_lumilevel.madx")

        print('Luminosities after leveling (crab cavities are not considered):')
        lumi.print_luminosity(mad, twiss_dfs,
                mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
                mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])
else:
    from scipy.optimize import least_squares

    print('Luminosities before leveling (crab cavities are not considered):')
    lumi.print_luminosity(mad, twiss_dfs,
            mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
            mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])

    # Leveling in IP8
    L_target_ip8 = mask_parameters['par_lumi_ip8']
    def function_to_minimize_ip8(sep8v_m):
        my_dict_IP8=lumi.get_luminosity_dict(
            mad, twiss_dfs, 'ip8', mask_parameters['par_nco_IP8'])
        my_dict_IP8['y_1']=np.abs(sep8v_m)
        my_dict_IP8['y_2']=-np.abs(sep8v_m)
        return np.abs(lumi.L(**my_dict_IP8) - L_target_ip8)
    sigma_x_b1_ip8=np.sqrt(twiss_dfs['lhcb1'].loc['ip8:1'].betx*mad.sequence.lhcb1.beam.ex)
    optres_ip8=least_squares(function_to_minimize_ip8, sigma_x_b1_ip8)
    mad.globals['on_sep8'] = np.sign(mad.globals['on_sep8']) * np.abs(optres_ip8['x'][0])*1e3

    # Halo collision in IP2
    sigma_y_b1_ip2=np.sqrt(twiss_dfs['lhcb1'].loc['ip2:1'].bety*mad.sequence.lhcb1.beam.ey)
    mad.globals['on_sep2']=np.sign(mad.globals['on_sep2'])*mask_parameters['par_fullsep_in_sigmas_ip2']*sigma_y_b1_ip2/2*1e3

    # Re-save knobs
    mad.input('exec, crossing_save')

    print('Luminosities after leveling (crab cavities are not considered):')
    lumi.print_luminosity(mad, twiss_dfs,
            mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
            mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])

if force_leveling is not None:
    for kk in force_leveling.keys():
        mad.globals[kk] = force_leveling[kk]
    mad.input('exec, crossing_save')

# Check machine after leveling
mad.input('exec, crossing_restore')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_after_leveling',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)



#####################
# Force on_disp = 0 #
#####################

mad.globals.on_disp = 0.
# will be restored later


###################################
# Compute beam-beam configuration #
###################################

# Prepare bb dataframes
if enable_bb_python:
    bb_dfs = pm.generate_bb_dataframes(mad,
        ip_names=['ip1', 'ip2', 'ip5', 'ip8'],
        harmonic_number=35640,
        numberOfLRPerIRSide=python_parameters['numberOfLRPerIRSide'],
        bunch_spacing_buckets=python_parameters['bunch_spacing_buckets'],
        numberOfHOSlices=python_parameters['numberOfHOSlices'],
        bunch_population_ppb=python_parameters['bunch_population_ppb'],
        sigmaz_m=python_parameters['sigmaz_m'],
        z_crab_twiss=python_parameters['z_crab_twiss'],
        remove_dummy_lenses=True)

    # Here the datafremes can be edited, e.g. to set bbb intensity


###################
# Generate beam 4 #
###################

if generate_b4_from_b2:
    mad_b4 = Madx()
    ost.build_sequence(mad_b4, beam=4)
    ost.apply_optics(mad_b4, optics_file=optics_file)

    pm.configure_b4_from_b2(mad_b4, mad)

    twiss_dfs_b2, other_data_b2 = ost.twiss_and_check(mad,
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


##################################################
# Select mad instance for tracking configuration #
##################################################

# We working exclusively on the sequence to track
# Select mad object
if track_from_b4_mad_instance:
    mad_track = mad_b4
else:
    mad_track = mad

mad_collider = mad
del(mad)

# Twiss machine to track
twiss_dfs, other_data = ost.twiss_and_check(mad_track, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_track_intermediate',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)


#####################
# Install bb lenses #
#####################

# Python approach
if enable_bb_python:
    if track_from_b4_mad_instance:
        bb_df_track = bb_dfs['b4']
        assert(sequence_to_track=='lhcb2')
    else:
        bb_df_track = bb_dfs['b1']
        assert(sequence_to_track=='lhcb1')

    pm.install_lenses_in_sequence(mad_track, bb_df_track, sequence_to_track)

    # Disable bb (to be activated later)
    mad_track.globals.on_bb_charge = 0
else:
    bb_df_track = None

# Legacy bb macros
if enable_bb_legacy:
    assert(beam_to_configure == 1)
    assert(not(track_from_b4_mad_instance))
    assert(not(enable_bb_python))
    mad_track.call("modules/module_03_beambeam.madx")


#########################
# Install crab cavities #
#########################
mad_track.call("optics_toolkit/enable_crabcavities.madx")
# They are left off, they will be swiched on at the end:
mad_track.globals.on_crab1 = 0
mad_track.globals.on_crab5 = 0

##############################################
# Save references for tuning and corrections #
##############################################
mad_track.call('modules/submodule_04_1b_save_references.madx')


#####################
# Force on_disp = 0 #
#####################

mad_track.globals.on_disp = 0.
# will be restored later


#############
# Final use #
#############

mad_track.use(sequence_to_track)
# Disable use
mad_track._use = mad_track.use
mad_track.use = None


##############################
# Install and correct errors #
##############################

if python_parameters['enable_multipolar_errors']:
    mad_track.call('modules/module_04_errors.madx')
else:
    # Synthesize knobs
    mad_track.call('modules/submodule_04a_s1_prepare_nom_twiss_table.madx')
    mad_track.input('exec, crossing_disable;')
    mad_track.call('modules/submodule_04e_s1_synthesize_knobs.madx')
    mad_track.input('exec, crossing_restore;')

###############################
# Machine tuning (enables bb) #
###############################

mad_track.call("modules/submodule_05a_MO.madx")

cmrskew_0 = mad_track.globals.cmrskew
cmiskew_0 = mad_track.globals.cmiskew

# Test old approach
mad_track.call('modules/submodule_05b_coupling.madx')
cmrskew_legacy = mad_track.globals.cmrskew
cmiskew_legacy = mad_track.globals.cmiskew
cta_legacy = pm.coupling_measurement(mad_track,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        sequence_name='lhcb1', skip_use=True)



# Test new approach
mad_track.globals.cmrskew = cmrskew_0
mad_track.globals.cmiskew = cmiskew_0

pm.coupling_correction(mad_track, n_iterations=2,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        cmr_knob_name = 'cmrskew', cmi_knob_name = 'cmiskew',
        sequence_name='lhcb1', skip_use=True)
cmrskew_new = mad_track.globals.cmrskew
cmiskew_new = mad_track.globals.cmiskew
cta_new = pm.coupling_measurement(mad_track,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        sequence_name='lhcb1', skip_use=True)

print(f'cmrskew_legacy = {cmrskew_legacy}')
print(f'cmrskew_new = {cmrskew_new}')
print(f'cmiskew_legacy = {cmiskew_legacy}')
print(f'cmiskew_new = {cmiskew_new}')
print(f'cta_legacy = {cta_legacy}')
print(f'cta_new = {cta_new}')

prrrrrr



# Correct linear coupling
qx_fractional, qx_integer = np.modf(mask_parameters['par_qx0'])
qy_fractional, qy_integer = np.modf(mask_parameters['par_qy0'])

mad_track.call("modules/submodule_05b_coupling.madx")
# DEBUG
qx_fractional = .31
qy_fractional = .32
beam_name = sequence_to_track[-2:]
pm.coupling_correction(mad_track, n_iterations=2,
        qx_integer=qx_integer, qy_integer=qy_integer,
        qx_fractional=qx_fractional, qy_fractional=qy_fractional,
        tune_knob1_name='kqtf.'+beam_name,
        tune_knob2_name='kqtd.'+beam_name,
        cmr_knob_name = 'cmrskew', cmi_knob_name = 'cmiskew',
        sequence_name=sequence_to_track, skip_use=True)


mad_track.call("modules/submodule_05c_limit.madx")
mad_track.call("modules/submodule_05d_matching.madx")
mad_track.call("modules/submodule_05e_corrvalue.madx")
mad_track.call("modules/submodule_05f_final.madx")

# Switch on crab cavities
mad_track.globals.on_crab1 = mad_track.globals.par_crab1
mad_track.globals.on_crab5 = mad_track.globals.par_crab5


#####################
# Generate sixtrack #
#####################

if enable_bb_legacy:
    mad_track.call("modules/module_06_generate.madx")
else:
    pm.generate_sixtrack_input(mad_track,
        seq_name=sequence_to_track,
        bb_df=bb_df_track,
        output_folder='./',
        reference_bunch_charge_sixtrack_ppb=(
            mad_track.sequence[sequence_to_track].beam.npart),
        emitnx_sixtrack_um=(
            mad_track.sequence[sequence_to_track].beam.exn),
        emitny_sixtrack_um=(
            mad_track.sequence[sequence_to_track].beam.eyn),
        sigz_sixtrack_m=(
            mad_track.sequence[sequence_to_track].beam.sigt),
        sige_sixtrack=(
            mad_track.sequence[sequence_to_track].beam.sige),
        ibeco_sixtrack=1,
        ibtyp_sixtrack=0,
        lhc_sixtrack=2,
        ibbc_sixtrack=0,
        radius_sixtrack_multip_conversion_mad=0.017,
        skip_mad_use=True)


#######################################
# Save optics and orbit at start ring #
#######################################

optics_orbit_start_ring = pm.get_optics_and_orbit_at_start_ring(
        mad_track, sequence_to_track, skip_mad_use=True)
with open('./optics_orbit_at_start_ring.pkl', 'wb') as fid:
    pickle.dump(optics_orbit_start_ring, fid)


#############################
# Generate pysixtrack lines #
#############################

if enable_bb_legacy:
    print('Pysixtrack line is not generated with bb legacy macros')
else:
    pysix_fol_name = "./pysixtrack"
    dct_pysxt = pm.generate_pysixtrack_line_with_bb(mad_track,
        sequence_to_track, bb_df_track,
        closed_orbit_method='from_mad',
        pickle_lines_in_folder=pysix_fol_name,
        skip_mad_use=True)

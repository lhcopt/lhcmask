import os
import sys
import pickle

import numpy as np

from config import python_parameters, mask_parameters
from config import knob_settings, knob_names


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
enable_lumi_control = python_parameters['enable_lumi_control']
enable_multipolar_errors = python_parameters['enable_multipolar_errors']
enable_crabs = python_parameters['enable_crabs']

# Make links
for kk in links.keys():
    if os.path.exists(kk):
        os.remove(kk)
    os.symlink(os.path.abspath(links[kk]), kk)

# Create empty temp folder
os.system('rm -r temp')
os.system('mkdir temp')

# Execute customization script if present
os.system('bash customization.bash')

# Import pymask
sys.path.append('./modules')
import pymask as pm

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

if not(enable_crabs):
    knob_settings['par_crab1'] = 0.
    knob_settings['par_crab5'] = 0.

########################
# Build MAD-X instance #
########################

# Start mad
Madx = pm.Madxp
mad = Madx(command_log="mad_collider.log")

# Build sequence (alse creates link to optics_toolkit and calls it)
ost.build_sequence(mad, beam=beam_to_configure)

# Set twiss formats for MAD-X parts (macro from opt. toolkit)
mad.input('exec, twiss_opt;')

# Apply optics
ost.apply_optics(mad, optics_file=optics_file)

# Pass parameters to mad
mad.set_variables_from_dict(params=mask_parameters)

# Prepare auxiliary mad variables
#mad.input("call, file='modules/submodule_01a_preparation.madx';")


# Attach beams to sequences
#mad.input("call, file='modules/submodule_01b_beam.madx';")

# Set energy
mad.globals.nrj = mask_parameters['par_beam_energy_tot']
for ss in mad.sequence.keys():
    if ss == 'lhcb1':
        beam_bv = 1
        bv_aux = 1
    elif ss == 'lhcb2':
        if int(beam_to_configure) == 4:
            ss_beam_bv = 1
            ss_bv_aux = -1
        else:
            ss_beam_bv = -1
            ss_bv_aux = 1
    mad.globals['bv_aux'] = ss_bv_aux
    mad.input(f'''
    beam, particle=proton,sequence={ff},
    energy={mask_parameters['par_beam_energy_tot']},
    sigt={mask_parameters['par_beam_sigt']},
    bv={ss_beam_bv},
    npart={mask_parameters['par_beam_npart'],
    sige={mask_parameters['par_beam_sige'],
    ex=epsx,ey=epsy;


# Test machine before any change
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_from_optics',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips,
        check_separations_at_ips=check_separations_at_ips)

# Set IP1-IP5 phase and store corresponding reference
mad.input("call, file='modules/submodule_01c_phase.madx';")

# Set optics-specific knobs
ost.set_optics_specific_knobs(mad, knob_settings, mode)

# Crossing-save and some reference measurements
mad.input('exec, crossing_save;')
mad.input("call, file='modules/submodule_01e_final.madx';")


#################################
# Check bahavior of orbit knobs #
#################################

# Check flat machine
mad.input('exec, crossing_disable;')
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
mad.input('exec, crossing_restore;')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_with_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)


#################################
# Set luminosity in IP2 and IP8 #
#################################

if len(sequences_to_check) == 2:
    print('Luminosities before leveling (crab cavities are not considered):')
    pm.print_luminosity(mad, twiss_dfs,
            mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
            mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])
else:
    print('Warning: Luminosity computation requires two beams')


if not enable_lumi_control:
    print('Separations in IP2 and IP8 are left untouched')
elif enable_bb_legacy or mode=='b4_without_bb':
    mad.use(f'lhcb{beam_to_configure}')
    if mode=='b4_without_bb':
        print('Leveling not working in this mode!')
    else:
        # Luminosity levelling
        mad.input("call, file='modules/module_02_lumilevel.madx';")
else:
    print('Start pythonic leveling:')
    ost.lumi_control(mad, twiss_dfs, python_parameters,
            mask_parameters, knob_names)

# Force leveling
if force_leveling is not None:
    for kk in force_leveling.keys():
        mad.globals[kk] = force_leveling[kk]

# Re-save knobs (for the last time!)
mad.input('exec, crossing_save;')

# Check machine after leveling
mad.input('exec, crossing_restore;')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_after_leveling',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips,
        check_separations_at_ips=check_separations_at_ips)

if len(sequences_to_check) == 2:
    print('Luminosities after leveling (crab cavities are not considered):')
    pm.print_luminosity(mad, twiss_dfs,
            mask_parameters['par_nco_IP1'], mask_parameters['par_nco_IP2'],
            mask_parameters['par_nco_IP5'], mask_parameters['par_nco_IP8'])
else:
    print('Luminosity computation requires two beams')


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
        z_crab_twiss=python_parameters['z_crab_twiss']*float(enable_crabs),
        remove_dummy_lenses=True)

    # Here the datafremes can be edited, e.g. to set bbb intensity


###################
# Generate beam 4 #
###################

if generate_b4_from_b2:
    mad_b4 = Madx(command_log="mad_b4.log")
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

# We will be working exclusively on the sequence to track
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
    mad_track.input("call, file='modules/module_03_beambeam.madx';")


#########################
# Install crab cavities #
#########################
if enable_crabs:
    mad_track.input("call, file='optics_toolkit/enable_crabcavities.madx';")
    # They are left off, they will be swiched on at the end:
    mad_track.globals.on_crab1 = 0
    mad_track.globals.on_crab5 = 0

##############################################
# Save references for tuning and corrections #
##############################################
mad_track.input("call, file='modules/submodule_04_1b_save_references.madx';")


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

if enable_multipolar_errors:
    mad_track.input("call, file='modules/module_04_errors.madx';")
else:
    # Synthesize knobs
    mad_track.input('call, file="modules/submodule_04a_s1_prepare_nom_twiss_table.madx";')
    if python_parameters['enable_knob_synthesis']:
        mad_track.input('exec, crossing_disable;')
        mad_track.input("call, file='modules/submodule_04e_s1_synthesize_knobs.madx';")
        mad_track.input('exec, crossing_restore;')

##################
# Machine tuning #
##################

# Enable bb for matchings
if mask_parameters['par_match_with_bb'] == 1:
    mad_track.globals['on_bb_charge'] = 1
else:
    mad_track.globals['on_bb_charge'] = 0

# Switch on octupoles
mad_track.input("call, file='modules/submodule_05a_MO.madx';")

# Correct linear coupling
qx_fractional, qx_integer = np.modf(mask_parameters['par_qx0'])
qy_fractional, qy_integer = np.modf(mask_parameters['par_qy0'])
coupl_corr_info = pm.coupling_correction(mad_track,
        n_iterations=python_parameters['N_iter_coupling'],
        qx_integer=qx_integer, qy_integer=qy_integer,
        qx_fractional=qx_fractional, qy_fractional=qy_fractional,
        tune_knob1_name=knob_names['qknob_1'][sequence_to_track],
        tune_knob2_name=knob_names['qknob_2'][sequence_to_track],
        cmr_knob_name=knob_names['cmrknob'][sequence_to_track],
        cmi_knob_name=knob_names['cmiknob'][sequence_to_track],
        sequence_name=sequence_to_track, skip_use=True)

# Add custom values to coupling knobs
mad_track.globals[knob_names['cmrknob'][sequence_to_track]] += python_parameters['delta_cmr']
mad_track.globals[knob_names['cmiknob'][sequence_to_track]] += python_parameters['delta_cmi']

# Check strength limits
if enable_multipolar_errors:
    mad_track.input('call, file="errors/HL-LHC/corr_limit.madx";')

# Rematch the orbit at IPs
mad_track.input("call, file='tools/rematchCOIP.madx';")

# Rematch the CO in the arc for dispersion correction
if mad_track.globals.on_disp != 0:
    mad_track.input("call, file='tools/rematchCOarc.madx';")

# Match tunes and chromaticities
pm.match_tune_and_chromaticity(mad_track,
        q1=mask_parameters['par_qx0'],
        q2=mask_parameters['par_qy0'],
        dq1=mask_parameters['par_chromaticity_x'],
        dq2=mask_parameters['par_chromaticity_y'],
        tune_knob1_name=knob_names['qknob_1'][sequence_to_track],
        tune_knob2_name=knob_names['qknob_2'][sequence_to_track],
        chromaticity_knob1_name=knob_names['chromknob_1'][sequence_to_track],
        chromaticity_knob2_name=knob_names['chromknob_2'][sequence_to_track],
        sequence_name=sequence_to_track, skip_use=True)

# Check strength limits
if enable_multipolar_errors:
    mad_track.input("call, file='errors/HL-LHC/corr_value_limit.madx';")

# Switch on bb lenses
mad_track.globals.on_bb_charge = 1.

# Switch on RF cavities
mad_track.globals['vrf400'] = mask_parameters['par_vrf_total']
if sequence_to_track == 'lhcb1':
    mad_track.globals['lagrf400.b1'] = 0.5
elif sequence_to_track == 'lhcb2':
    mad_track.globals['lagrf400.b2'] = 0.

# Switch on crab cavities
if enable_crabs:
    mad_track.globals.on_crab1 = knob_settings['on_crab1']
    mad_track.globals.on_crab5 = knob_settings['on_crab5']

#####################
# Generate sixtrack #
#####################

if enable_bb_legacy:
    mad_track.input("call, file='modules/module_06_generate.madx'")
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

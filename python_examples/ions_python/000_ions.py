import os
import sys
import pickle

import numpy as np


#####################################################
# Read general configurations and setup envirnoment #
#####################################################

from config import configuration

mode = configuration['mode']
tol_beta = configuration['tol_beta']
tol_sep = configuration['tol_sep']
flat_tol = configuration['tol_co_flatness']
links = configuration['links']
optics_file = configuration['optics_file']
check_betas_at_ips = configuration['check_betas_at_ips']
check_separations_at_ips = configuration['check_separations_at_ips']
save_intermediate_twiss = configuration['save_intermediate_twiss']
enable_lumi_control = configuration['enable_lumi_control']
enable_imperfections = configuration['enable_imperfections']
enable_crabs = configuration['enable_crabs']
match_q_dq_with_bb = configuration['match_q_dq_with_bb']
knob_settings = configuration['knob_settings']
knob_names = configuration['knob_names']


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

# Define configuration
(beam_to_configure, sequences_to_check, sequence_to_track, generate_b4_from_b2,
    track_from_b4_mad_instance, enable_bb_python, enable_bb_legacy,
    force_disable_check_separations_at_ips,
    ) = pm.get_pymask_configuration(mode)

if force_disable_check_separations_at_ips:
    check_separations_at_ips = False

if not(enable_crabs):
    knob_settings['par_crab1'] = 0.
    knob_settings['par_crab5'] = 0.

########################
# Build MAD-X instance #
########################

# Start mad
Madx = pm.Madxp
mad = Madx(command_log="mad_collider.log")

# Set verbose flag
mad.globals.par_verbose = int(configuration['verbose_mad_parts'])

# Build sequence (alse creates link to optics_toolkit and calls it)
ost.build_sequence(mad, beam=beam_to_configure)

# Set twiss formats for MAD-X parts (macro from opt. toolkit)
mad.input('exec, twiss_opt;')

# Apply optics
ost.apply_optics(mad, optics_file=optics_file)

# Attach beam to sequences
mad.globals.nrj = configuration['beam_energy_tot']
gamma_rel = configuration['beam_energy_tot']/configuration['mass']
for ss in mad.sequence.keys():
    # bv and bv_aux flags
    if ss == 'lhcb1':
        ss_beam_bv, ss_bv_aux = 1, 1
    elif ss == 'lhcb2':
        if int(beam_to_configure) == 4:
            ss_beam_bv, ss_bv_aux = 1, -1
        else:
            ss_beam_bv, ss_bv_aux = -1, 1

    mad.globals['bv_aux'] = ss_bv_aux
    mad.input(f'''
    beam, particle=ion,sequence={ss},
        energy={configuration['beam_energy_tot']},
        sigt={configuration['beam_sigt']},
        bv={ss_beam_bv},
        npart={configuration['beam_npart']},
        sige={configuration['beam_sige']},
        ex={configuration['beam_norm_emit_x'] * 1e-6 / gamma_rel},
        ey={configuration['beam_norm_emit_y'] * 1e-6 / gamma_rel},
        mass = {configuration['mass']},
        charge={configuration['charge']},
    ''')


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
            configuration['nco_IP1'], configuration['nco_IP2'],
            configuration['nco_IP5'], configuration['nco_IP8'])
else:
    print('Warning: Luminosity computation requires two beams')


if not enable_lumi_control:
    print('Separations in IP2 and IP8 are left untouched')
elif enable_bb_legacy or mode=='b4_without_bb':
    mad.use(f'lhcb{beam_to_configure}')
    if mode=='b4_without_bb':
        print('Leveling not working in this mode!')
    else:
        # Modified module 2 for offset levelling in IP1,2,5,8
        vars_for_legacy_level = ['lumi_ip8',
            'nco_IP1', 'nco_IP2', 'nco_IP5', 'nco_IP8']
        mad.set_variables_from_dict({
            'par_'+kk: configuration[kk] for kk in vars_for_legacy_level})
        mad.input('''
        print, text="";
        print, text="";
        print, text="+++++++++++++++++++++++++++++++";
        print, text="++ START MODULE 2: LEVELLING ++";
        print, text="+++++++++++++++++++++++++++++++";
        print, text="";
        print, text="";

        call, file="beambeam_macros/macro_bb.madx";                  ! macros for beam-beam

        exec, DEFINE_BB_PARAM;  !Define main beam-beam parameters

        !Switch on Xscheme in precollision
        on_disp:=0;
        halo1=0;halo2=0;halo5=0;halo8=0;  !halo collision at 5 sigma's in Alice
        !number of collision/turn at IP1/2/5/8 
        nco_IP1 = par_nco_IP1;
        nco_IP5 = par_nco_IP5;
        nco_IP2 = par_nco_IP2;
        nco_IP8 = par_nco_IP8;
        exec, LEVEL_PARALLEL_OFFSET_FOR(par_lumi_ip8, 8); value,halo8;
        exec, LEVEL_PARALLEL_OFFSET_FOR(par_lumi_ip1, 1); value,halo1;
        exec, LEVEL_PARALLEL_OFFSET_FOR(par_lumi_ip2, 2); value,halo2;
        exec, LEVEL_PARALLEL_OFFSET_FOR(par_lumi_ip5, 5); value,halo5;
        !Redefine the on_sep's accordingly
        exec, CALCULATE_XSCHEME(halo1,halo2,halo5,halo8);
        ! Saving new crossing scheme with separation
        on_disp=on_dispaux; ! reset on_disp before saving
        exec, crossing_save;

        if (mylhcbeam==1) { use, sequence=lhcb1; } else { use, sequence=lhcb2; };

        ''')        
else:
    print('Start pythonic leveling:')
    ost.lumi_control(mad, twiss_dfs, configuration, knob_names)

# Force leveling
if 'force_leveling' in configuration.keys():
    force_leveling  = configuration['force_leveling']
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
            configuration['nco_IP1'], configuration['nco_IP2'],
            configuration['nco_IP5'], configuration['nco_IP8'])
else:
    print('Luminosity computation requires two beams')


#####################
# Force on_disp = 0 #
#####################

mad.globals.on_disp = 0.
# will be restored later


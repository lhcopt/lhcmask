import os
import json
import yaml
import numpy as np
import pymask as pm
import xobjects as xo
import xtrack as xt

from pymask.line_preparation import rename_coupling_knobs_and_coefficients
from pymask.line_preparation import define_octupole_current_knobs
from pymask.line_preparation import add_correction_term_to_dipole_correctors

# Import user-defined optics-specific tools
import optics_specific_tools as ost

# Read config file
with open('config.yaml','r') as fid:
    configuration = yaml.safe_load(fid)

# Make links
if configuration['links']['tracking_tools'] == 'auto':
    configuration['links']['tracking_tools'] = str(pm._pkg_root.parent.parent.absolute())

for kk in configuration['links'].keys():
    os.system(f'rm {kk}')
    os.symlink(os.path.abspath(configuration['links'][kk]), kk)

# Create empty temp folder
os.system('rm -r temp; mkdir temp')

# Set crab knobs to zero if crabbing is disabled
if not(configuration['enable_crabs']):
    configuration['knob_settings']['par_crab1'] = 0.
    configuration['knob_settings']['par_crab5'] = 0.

# Start mad
Madx = pm.Madxp
mad = Madx(command_log="mad_collider.log")
mad.globals.par_verbose = int(configuration['verbose_mad_parts'])

# Build sequence (also creates link to optics_toolkit and calls it)
ost.build_sequence(mad, beam=1, configuration=configuration)

# Set twiss formats for MAD-X parts (macro from opt. toolkit)
mad.input('exec, twiss_opt;')

# Apply optics
ost.apply_optics(mad, optics_file=configuration['optics_file'])

# Attach beam to sequences
pm.attach_beam_to_sequences(mad, beam_configuration=configuration)

mad.use('lhcb1')
mad.twiss()
mad.use('lhcb2')
mad.twiss()

# Generate beam 4
mad_b4 = Madx(command_log="mad_b4.log")
ost.build_sequence(mad_b4, beam=4, configuration=configuration)
pm.configure_b4_from_b2(mad_b4, mad)

lines_co_ref = {}
lines_co_ref['lhcb1_co_ref'] = xt.Line.from_madx_sequence(mad.sequence.lhcb1,
    deferred_expressions=True,
    expressions_for_element_types=('kicker', 'hkicker', 'vkicker'))
lines_co_ref['lhcb2_co_ref'] = xt.Line.from_madx_sequence(mad_b4.sequence.lhcb2,
    deferred_expressions=True,
    expressions_for_element_types=('kicker', 'hkicker', 'vkicker'))


# Set IP1-IP5 phase and store corresponding reference
mad.input("call, file='modules/submodule_01c_phase.madx';")

# Set optics-specific knobs
ost.set_optics_specific_knobs(mad, configuration['knob_settings'])

# Crossing-save and some reference measurements
mad.input('exec, crossing_save;')
mad.input("call, file='modules/submodule_01e_final.madx';")

# Check flat machine
mad.input('exec, crossing_disable;')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check=['lhcb1', 'lhcb2'],
        tol_beta=configuration['tol_beta'], tol_sep=configuration['tol_sep'],
        twiss_fname='twiss_no_crossing', save_twiss_files=True,
        check_betas_at_ips=configuration['check_betas_at_ips'],
        check_separations_at_ips=configuration['check_separations_at_ips'])

# Check orbit flatness
for ss in twiss_dfs.keys():
    tt = twiss_dfs[ss]
    assert np.max(np.abs(tt.x)) < configuration['tol_co_flatness']
    assert np.max(np.abs(tt.y)) < configuration['tol_co_flatness']

# Check machine after crossing restore
mad.input('exec, crossing_restore;')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check=['lhcb1', 'lhcb2'],
        tol_beta=configuration['tol_beta'], tol_sep=configuration['tol_sep'],
        twiss_fname='twiss_with_crossing', save_twiss_files=True,
        check_betas_at_ips=configuration['check_betas_at_ips'],
        check_separations_at_ips=configuration['check_separations_at_ips'])

# Set luminosity in IP2 and IP8
if configuration['enable_lumi_control']:
    ost.lumi_control(mad, twiss_dfs, configuration, configuration['knob_names'])

# Re-save orbit knobs knobs
mad.input('exec, crossing_save;')

# Force on_disp = 0
mad.globals.on_disp = 0. # will be restored later

# Update beam 4
pm.configure_b4_from_b2(mad_b4, mad)

lines_to_track = {}

for sequence_to_track, mad_track in zip(['lhcb1', 'lhcb2'], [mad, mad_b4]):

    # Install crab cavities
    if configuration['enable_crabs']:
        mad_track.input("call, file='optics_toolkit/enable_crabcavities.madx';")
        # They are left off, they will be swiched on at the end
        mad_track.globals.on_crab1 = 0
        mad_track.globals.on_crab5 = 0

    # Save references for tuning and corrections (does crossing restore, restores on_disp)
    mad_track.input("call, file='modules/submodule_04_1b_save_references.madx';")

    # Force on_disp = 0
    mad_track.globals.on_disp = 0. # will be restored later

    # Final use --> disable use
    mad_track.use(sequence_to_track)
    mad_track._use = mad_track.use
    mad_track.use = None

    # Install and correct errors
    if configuration['enable_imperfections']:
        mad_track.set_variables_from_dict(
                configuration['pars_for_imperfections'])
        mad_track.input("call, file='modules/module_04_errors.madx';")
    else:
        # Synthesize knobs
        mad_track.input('call, file="modules/submodule_04a_s1_prepare_nom_twiss_table.madx";')
        if configuration['enable_knob_synthesis']:
            mad_track.input('exec, crossing_disable;')
            mad_track.input("call, file='modules/submodule_04e_s1_synthesize_knobs.madx';")
        mad_track.input('exec, crossing_restore;')

    # Rename coupling knobs

    # Switch on octupoles
    brho = mad_track.globals.nrj*1e9/mad_track.globals.clight
    i_oct = configuration['oct_current']
    beam_str = {'lhcb1':'b1', 'lhcb2':'b2'}[sequence_to_track]
    for ss in '12 23 34 45 56 67 78 81'.split():
        mad_track.input(f'kof.a{ss}{beam_str} = kmax_mo*({i_oct})/imax_mo/({brho});')
        mad_track.input(f'kod.a{ss}{beam_str} = kmax_mo*({i_oct})/imax_mo/({brho});')

    # Correct linear coupling
    qx_fractional, qx_integer = np.modf(configuration['qx0'])
    qy_fractional, qy_integer = np.modf(configuration['qy0'])
    knob_names = configuration['knob_names']
    mad_track.globals[knob_names['cmrknob'][sequence_to_track]] = 0.
    mad_track.globals[knob_names['cmiknob'][sequence_to_track]] = 0.
    coupl_corr_info = pm.coupling_correction(mad_track,
            n_iterations=configuration['N_iter_coupling'],
            qx_integer=qx_integer, qy_integer=qy_integer,
            qx_fractional=qx_fractional, qy_fractional=qy_fractional,
            tune_knob1_name=knob_names['qknob_1'][sequence_to_track],
            tune_knob2_name=knob_names['qknob_2'][sequence_to_track],
            cmr_knob_name=knob_names['cmrknob'][sequence_to_track],
            cmi_knob_name=knob_names['cmiknob'][sequence_to_track],
            sequence_name=sequence_to_track, skip_use=True)

    # Add custom values to coupling knobs
    mad_track.globals[knob_names['cmrknob'][sequence_to_track]] += configuration['delta_cmr']
    mad_track.globals[knob_names['cmiknob'][sequence_to_track]] += configuration['delta_cmi']

    # # Check strength limits
    # if enable_imperfections:
    #     mad_track.input('call, file="errors/HL-LHC/corr_limit.madx";')

    # Rematch the orbit at IPs and in the arcs
    mad_track.input("call, file='tools/rematchCOIP.madx';")
    if mad_track.globals.on_disp != 0: # to clean leakage from dispersion correction
        mad_track.input("call, file='tools/rematchCOarc.madx';")

    # Match tunes and chromaticities
    pm.match_tune_and_chromaticity(mad_track,
        q1=configuration['qx0'], q2=configuration['qy0'],
        dq1=configuration['chromaticity_x'], dq2=configuration['chromaticity_y'],
        tune_knob1_name=configuration['knob_names']['qknob_1'][sequence_to_track],
        tune_knob2_name=configuration['knob_names']['qknob_2'][sequence_to_track],
        chromaticity_knob1_name=configuration['knob_names']['chromknob_1'][sequence_to_track],
        chromaticity_knob2_name=configuration['knob_names']['chromknob_2'][sequence_to_track],
        sequence_name=sequence_to_track, skip_use=True)

    # Switch on RF cavities
    mad_track.globals['vrf400'] = configuration['vrf_total']
    mad_track.globals['lagrf400.b1'] = 0.5
    mad_track.globals['lagrf400.b2'] = 0.

    # Switch on crab cavities
    if configuration['enable_crabs']:
        mad_track.globals.on_crab1 = configuration['knob_settings']['on_crab1']
        mad_track.globals.on_crab5 = configuration['knob_settings']['on_crab5']

    # Generate xtrack line
    tracker, line_bb_for_tracking_dict = pm.generate_xsuite_line(
                        mad_track, sequence_to_track)

    # Some fixes on knob definitions
    line = tracker.line
    rename_coupling_knobs_and_coefficients(line=line,
                                           beamn=int(sequence_to_track[-1]))
    define_octupole_current_knobs(line=line, beamn=int(sequence_to_track[-1]))
    lines_to_track[sequence_to_track] = line
    # Set octupoles
    line.vars['i_oct_b1'] = i_oct
    line.vars['i_oct_b2'] = i_oct

    # Save xtrack line to json
    with open('xsuite_line_' + sequence_to_track.replace('b2', 'b4') + '.json', 'w') as fid:
        json.dump(line_bb_for_tracking_dict, fid, cls=xo.JEncoder)

    # Save mad sequence
    pm.save_mad_sequence_and_error(mad_track, sequence_to_track,
        filename=sequence_to_track.replace('b2', 'b4'))

collider = xt.Multiline(
    lines={
        'lhcb1': lines_to_track['lhcb1'],
        'lhcb2': lines_to_track['lhcb2'],
        'lhcb1_co_ref': lines_co_ref['lhcb1_co_ref'],
        'lhcb2_co_ref': lines_co_ref['lhcb2_co_ref'],
    })

collider['lhcb1_co_ref'].particle_ref = collider['lhcb1'].particle_ref.copy()
collider['lhcb2_co_ref'].particle_ref = collider['lhcb2'].particle_ref.copy()

# Save the two lines to json
with open('collider.json', 'w') as fid:
    dct = collider.to_dict()
    json.dump(dct, fid, cls=xo.JEncoder)

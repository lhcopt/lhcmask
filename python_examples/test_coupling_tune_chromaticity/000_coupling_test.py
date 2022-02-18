import os
import sys
import json
import yaml
import numpy as np

#####################################################
# Read general configurations and setup envirnoment #
#####################################################

assert not(os.path.isfile('config.yaml')
           and os.path.isfile('config.py')), (
    "Please specify only a config file (yaml or py)")

try:
    with open('config.yaml','r') as fid:
        configuration = yaml.safe_load(fid)
except:
    from config import configuration

# Start tree_maker logging if log_file is present in config
try:
    import tree_maker
    if 'log_file' not in configuration.keys():
        tree_maker=None
except:
    tree_maker=None

if tree_maker is not None:
    tree_maker.tag_json.tag_it(configuration['log_file'], 'started')

mode = configuration['mode']
tol_beta = configuration['tol_beta']
tol_sep = configuration['tol_sep']
flat_tol = configuration['tol_co_flatness']
links = configuration['links']
optics_file = configuration['optics_file']
check_betas_at_ips = configuration['check_betas_at_ips']
check_separations_at_ips = configuration['check_separations_at_ips']
save_intermediate_twiss = configuration['save_intermediate_twiss']
knob_settings = configuration['knob_settings']

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
import pymask.coupling as pc

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

# Attach beam to sequences
mad.globals.nrj = configuration['beam_energy_tot']
gamma_rel = configuration['beam_energy_tot']/mad.globals.pmass
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
    beam, particle=proton,sequence={ss},
        energy={configuration['beam_energy_tot']},
        sigt={configuration['beam_sigt']},
        bv={ss_beam_bv},
        npart={configuration['beam_npart']},
        sige={configuration['beam_sige']},
        ex={configuration['beam_norm_emit_x'] * 1e-6 / gamma_rel},
        ey={configuration['beam_norm_emit_y'] * 1e-6 / gamma_rel};
    ''')

# Set optics-specific knobs
ost.set_optics_specific_knobs(mad, knob_settings, mode)

# Synthesisze knobs
mad.call('modules/submodule_04_1b_save_references.madx')
mad.call('modules/submodule_04a_s1_prepare_nom_twiss_table.madx')
mad.call('modules/submodule_04e_s1_synthesize_knobs.madx')


cmrskew_test = 1e-4
cmiskew_test = 0.

# Introduce large coupling for testing
mad.globals.cmrskew = cmrskew_test
mad.globals.cmiskew = cmiskew_test

# Test old approach
mad.globals.qx0 = 62.31
mad.globals.qy0 = 60.32
mad.globals.qx00 = 62.
mad.globals.qy00 = 60.
mad.call('modules/submodule_05b_coupling.madx')
cmrskew_legacy = mad.globals.cmrskew
cmiskew_legacy = mad.globals.cmiskew
cta_legacy = pc.coupling_measurement(mad,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        sequence_name='lhcb1', skip_use=False)

# Test new approach
mad.globals.cmrskew = cmrskew_test
mad.globals.cmiskew = cmiskew_test

pc.coupling_correction(mad, n_iterations=2,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        cmr_knob_name = 'cmrskew', cmi_knob_name = 'cmiskew',
        sequence_name='lhcb1', skip_use=False)
cmrskew_new = mad.globals.cmrskew
cmiskew_new = mad.globals.cmiskew
cta_new = pc.coupling_measurement(mad,
        qx_integer=62., qy_integer=60.,
        qx_fractional=.31, qy_fractional=.32,
        tune_knob1_name='kqtf.b1', tune_knob2_name='kqtd.b1',
        sequence_name='lhcb1', skip_use=False)

print(f'cmrskew_legacy = {cmrskew_legacy}')
print(f'cmrskew_new = {cmrskew_new}')
print(f'cmiskew_legacy = {cmiskew_legacy}')
print(f'cmiskew_new = {cmiskew_new}')
print(f'cta_legacy = {cta_legacy}')
print(f'cta_new = {cta_new}')


pm.match_tune_and_chromaticity(mad,
        q1=62.30,
        q2=60.25,
        dq1=12.,
        dq2=14.,
        tune_knob1_name='kqtf.b1',
        tune_knob2_name='kqtd.b1',
        chromaticity_knob1_name='ksf.b1',
        chromaticity_knob2_name='ksd.b1',
        sequence_name='lhcb1', skip_use=False)

mad.twiss()

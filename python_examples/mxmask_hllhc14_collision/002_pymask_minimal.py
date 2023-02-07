import os
import json
import yaml
import pymask as pm
import xobjects as xo
import xtrack as xt
import xpart as xp

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
pm.attach_beam_to_sequences(mad, configuration=configuration)

# Set cavity phases
mad.globals['lagrf400.b1'] = 0.5
mad.globals['lagrf400.b2'] = 0.

mad.use('lhcb1')
mad.twiss()
mad.use('lhcb2')
mad.twiss()

# Generate beam 4
mad_b4 = Madx(command_log="mad_b4.log")
ost.build_sequence(mad_b4, beam=4, configuration=configuration)
pm.configure_b4_from_b2(mad_b4, mad)

lines_co_ref = {}
# lines_co_ref['lhcb1_co_ref'] = xt.Line.from_madx_sequence(mad.sequence.lhcb1,
#     deferred_expressions=True,
#     expressions_for_element_types=('kicker', 'hkicker', 'vkicker'))
lines_co_ref['lhcb2_co_ref'] = xt.Line.from_madx_sequence(mad_b4.sequence.lhcb2,
    deferred_expressions=True,
    expressions_for_element_types=('kicker', 'hkicker', 'vkicker'))

# Set IP1-IP5 phase and store corresponding reference
#mad.input("call, file='modules/submodule_01c_phase.madx';")

# # Update beam 4
# pm.configure_b4_from_b2(mad_b4, mad)

lines_to_track = {}
for sequence_to_track, mad_track in zip(['lhcb1', 'lhcb2'], [mad, mad_b4]):

    # # Install crab cavities
    # if configuration['enable_crabs']:
    #     mad_track.input("call, file='optics_toolkit/enable_crabcavities.madx';")
    #     # They are left off, they will be swiched on at the end
    #     mad_track.globals.on_crab1 = 0
    #     mad_track.globals.on_crab5 = 0

    # # Save references for tuning and corrections (does crossing restore, restores on_disp)
    # mad_track.input("call, file='modules/submodule_04_1b_save_references.madx';")

    # # Force on_disp = 0
    # mad_track.globals.on_disp = 0. # will be restored later

    # # Final use --> disable use
    # mad_track.use(sequence_to_track)

    # # Install and correct errors
    # if configuration['enable_imperfections']:
    #     mad_track.set_variables_from_dict(
    #             configuration['pars_for_imperfections'])
    #     mad_track.input("call, file='modules/module_04_errors.madx';")
    # else:
    #     # Synthesize knobs
    #     mad_track.input('call, file="modules/submodule_04a_s1_prepare_nom_twiss_table.madx";')
    #     if configuration['enable_knob_synthesis']:
    #         mad_track.input('exec, crossing_disable;')
    #         mad_track.input("call, file='modules/submodule_04e_s1_synthesize_knobs.madx';")
    #     mad_track.input('exec, crossing_restore;')

    # Prepare xsuite line
    line = xt.Line.from_madx_sequence(
        mad_track.sequence[sequence_to_track], apply_madx_errors=True,
        deferred_expressions=True)
    mad_beam = mad_track.sequence[sequence_to_track].beam
    line.particle_ref = xp.Particles(
        p0c = mad_beam.pc*1e9,
        q0 = mad_beam.charge,
        mass0 = mad_beam.mass*1e9,
    )
    # rename_coupling_knobs_and_coefficients(line=line,
    #                                        beamn=int(sequence_to_track[-1]))
    # define_octupole_current_knobs(line=line, beamn=int(sequence_to_track[-1]))
    lines_to_track[sequence_to_track] = line


collider = xt.Multiline(
    lines={
        'lhcb1': lines_to_track['lhcb1'],
        'lhcb2': lines_to_track['lhcb2'],
        #'lhcb1_co_ref': lines_co_ref['lhcb1_co_ref'],
        'lhcb2_co_ref': lines_co_ref['lhcb2_co_ref'],
    })

#collider['lhcb1_co_ref'].particle_ref = collider['lhcb1'].particle_ref.copy()
collider['lhcb2_co_ref'].particle_ref = collider['lhcb2'].particle_ref.copy()



# add_correction_term_to_dipole_correctors(collider)

# Save the two lines to json
with open('collider.json', 'w') as fid:
    dct = collider.to_dict()
    json.dump(dct, fid, cls=xo.JEncoder)

#### Debug: find issue with errors in b4

collider.build_trackers()

p_test = collider['lhcb2_co_ref'].build_particles(x=1e-4, y=1e-4)

collider['lhcb2'].track(particles=p_test.copy(), turn_by_turn_monitor='ONE_TURN_EBE')
mon_test = collider['lhcb2'].record_last_track

collider['lhcb2_co_ref'].track(particles=p_test.copy(), turn_by_turn_monitor='ONE_TURN_EBE')
mon_ref = collider['lhcb2_co_ref'].record_last_track



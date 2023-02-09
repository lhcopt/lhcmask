import os
import json
import yaml
import pymask as pm
import xobjects as xo
import xtrack as xt
import xpart as xp

from pymask.line_preparation import make_mad_environment
from pymask.line_preparation import rename_coupling_knobs_and_coefficients
from pymask.line_preparation import define_octupole_current_knobs
from pymask.line_preparation import add_correction_term_to_dipole_correctors
from pymask.line_preparation import install_correct_errors_and_synthesisize_knobs
from pymask.line_preparation import  save_lines_for_closed_orbit_reference

# Import user-defined optics-specific tools
import optics_specific_tools as ost

# Read config file
with open('config.yaml','r') as fid:
    configuration = yaml.safe_load(fid)

# Make mad environment
make_mad_environment(links=configuration['links'])

# Start mad
Madx = pm.Madxp
mad = Madx(command_log="mad_collider.log")
mad.globals.par_verbose = int(configuration['verbose_mad_parts'])

# Build sequence, load optics, define beam
ost.build_sequence(mad, beam=1, optics_version=configuration['optics_version'])
ost.apply_optics(mad, optics_file=configuration['optics_file'])
pm.attach_beam_to_sequences(mad, beam_configuration=configuration)

# Warm up :-)
mad.use('lhcb1'); mad.twiss(); mad.use('lhcb2'); mad.twiss()

# Generate beam 4
mad_b4 = Madx(command_log="mad_b4.log")
ost.build_sequence(mad_b4, beam=4, optics_version=configuration['optics_version'])
pm.configure_b4_from_b2(mad_b4, mad)

# Save lines for closed orbit reference
lines_co_ref = save_lines_for_closed_orbit_reference(mad, mad_b4)

# Set IP1-IP5 phase and store corresponding reference
mad.input("call, file='modules/submodule_01c_phase.madx';")
pm.configure_b4_from_b2(mad_b4, mad) # Update b4

lines_to_track = {}
for sequence_to_track, mad_track in zip(['lhcb1', 'lhcb2'], [mad, mad_b4]):

    # Final use
    mad_track.use(sequence_to_track)

    # We work exclusively on the flat machine
    mad_track.input('exec, crossing_disable;')
    mad_track.input('exec, crossing_save;') # In this way crossing_restore keeps the flat machine

    # Install and correct errors
    install_correct_errors_and_synthesisize_knobs(mad_track,
        enable_imperfections=configuration['enable_imperfections'],
        enable_knob_synthesis= configuration['enable_knob_synthesis'],
        pars_for_imperfections=configuration['pars_for_imperfections'])

    # Prepare xsuite line
    line = xt.Line.from_madx_sequence(
        mad_track.sequence[sequence_to_track], apply_madx_errors=True,
        deferred_expressions=True,
        replace_in_expr={'bv_aux': 'bvaux_'+sequence_to_track})
    mad_beam = mad_track.sequence[sequence_to_track].beam
    line.particle_ref = xp.Particles(p0c = mad_beam.pc*1e9,
        q0 = mad_beam.charge, mass0 = mad_beam.mass*1e9)

    # Prepare coupling and octupole knobs
    rename_coupling_knobs_and_coefficients(line=line,
                                           beamn=int(sequence_to_track[-1]))
    define_octupole_current_knobs(line=line, beamn=int(sequence_to_track[-1]))
    lines_to_track[sequence_to_track] = line


collider = xt.Multiline(
    lines={
        'lhcb1': lines_to_track['lhcb1'],
        'lhcb2': lines_to_track['lhcb2'],
        'lhcb1_co_ref': lines_co_ref['lhcb1_co_ref'],
        'lhcb2_co_ref': lines_co_ref['lhcb2_co_ref'],
    })

collider['lhcb1_co_ref'].particle_ref = collider['lhcb1'].particle_ref.copy()
collider['lhcb2_co_ref'].particle_ref = collider['lhcb2'].particle_ref.copy()

add_correction_term_to_dipole_correctors(collider)

# Save the two lines to json
with open('collider.json', 'w') as fid:
    dct = collider.to_dict()
    json.dump(dct, fid, cls=xo.JEncoder)

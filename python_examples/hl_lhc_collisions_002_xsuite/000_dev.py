import json
import numpy as np
import xtrack as xt
import xpart as xp
import xfields as xf

from pymask.beambeam import generate_set_of_bb_encounters_1beam
from pymask.beambeam import get_geometry_and_optics_b1_b2
from pymask.beambeam import get_survey_ip_position_b1_b2
from pymask.beambeam import get_partner_corrected_position_and_optics
from pymask.beambeam import compute_separations
from pymask.beambeam import compute_dpx_dpy
from pymask.beambeam import compute_local_crossing_angle_and_plane
from pymask.beambeam import compute_xma_yma
from pymask.beambeam import get_counter_rotating
from pymask.beambeam import crabbing_strong_beam

with open('../hl_lhc_collisions_000_b1_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b1 = json.load(fid)
with open('../hl_lhc_collisions_001_b4_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b4 = json.load(fid)
line_b1 = xt.Line.from_dict(dct_b1)
line_b4 = xt.Line.from_dict(dct_b4)


line_b1.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)
line_b4.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)


ip_names=['ip1', 'ip2', 'ip5', 'ip8']
numberOfLRPerIRSide=[25, 20, 25, 20]
harmonic_number=35640
bunch_spacing_buckets=10
numberOfHOSlices=11
bunch_num_particles=2.2e11
sigmaz_m = 0.076
nemitt_x = 2.5e-6
nemitt_y = 2.5e-6

circumference = line_b1.get_length()


# TODO: use keyword arguments
# TODO: what happens if bunch length is different for the two beams
bb_df_b1 = generate_set_of_bb_encounters_1beam(
    circumference, harmonic_number,
    bunch_spacing_buckets,
    numberOfHOSlices,
    bunch_num_particles, line_b1.particle_ref.q0,
    sigmaz_m, line_b1.particle_ref.beta0[0], ip_names, numberOfLRPerIRSide,
    beam_name = 'b1',
    other_beam_name = 'b2')


bb_df_b2 = generate_set_of_bb_encounters_1beam(
    circumference, harmonic_number,
    bunch_spacing_buckets,
    numberOfHOSlices,
    bunch_num_particles, line_b4.particle_ref.q0,
    sigmaz_m,
    line_b4.particle_ref.beta0[0], ip_names, numberOfLRPerIRSide,
    beam_name = 'b2',
    other_beam_name = 'b1')
bb_df_b2['atPosition'] = -bb_df_b2['atPosition'] # I am installing in b4 not in b2

from temp_module import install_dummy_bb_lenses

install_dummy_bb_lenses(bb_df=bb_df_b1, line=line_b1)
install_dummy_bb_lenses(bb_df=bb_df_b2, line=line_b4)

tracker_b1 = line_b1.build_tracker()
tracker_b4 = line_b4.build_tracker()

twiss_b1 = tracker_b1.twiss()
twiss_b2 = tracker_b4.twiss(reverse=True)

survey_b1 = tracker_b1.survey()
survey_b2 = tracker_b4.survey(reverse=True)

sigmas_b1 = twiss_b1.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
sigmas_b2 = twiss_b2.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)

# Use survey and twiss to get geometry and locations of all encounters
get_geometry_and_optics_b1_b2(
    mad=None,
    bb_df_b1=bb_df_b1,
    bb_df_b2=bb_df_b2,
    xsuite_line_b1=line_b1,
    xsuite_line_b2=line_b4,
    xsuite_twiss_b1=twiss_b1,
    xsuite_twiss_b2=twiss_b2,
    xsuite_survey_b1=survey_b1,
    xsuite_survey_b2=survey_b2,
    xsuite_sigmas_b1=sigmas_b1,
    xsuite_sigmas_b2=sigmas_b2,
)

# Get the position of the IPs in the surveys of the two beams
ip_position_df = get_survey_ip_position_b1_b2(mad=None, ip_names=ip_names,
    xsuite_survey_b1=survey_b1, xsuite_survey_b2=survey_b2)

# Get geometry and optics at the partner encounter
get_partner_corrected_position_and_optics(
        bb_df_b1, bb_df_b2, ip_position_df)

# Compute separation, crossing plane rotation, crossing angle and xma
for bb_df in [bb_df_b1, bb_df_b2]:
    compute_separations(bb_df)
    compute_dpx_dpy(bb_df)
    compute_local_crossing_angle_and_plane(bb_df)
    compute_xma_yma(bb_df)

# Get bb dataframe and mad model (with dummy bb) for beam 3 and 4
bb_df_b3 = get_counter_rotating(bb_df_b1)
bb_df_b4 = get_counter_rotating(bb_df_b2)
#generate_mad_bb_info(bb_df_b3, mode='dummy')
#generate_mad_bb_info(bb_df_b4, mode='dummy')

# Generate mad info
# generate_mad_bb_info(bb_df_b1, mode='from_dataframe',
#         madx_reference_bunch_num_particles=madx_reference_bunch_num_particles)
# generate_mad_bb_info(bb_df_b2, mode='from_dataframe',
#         madx_reference_bunch_num_particles=madx_reference_bunch_num_particles)
# generate_mad_bb_info(bb_df_b3, mode='from_dataframe',
#         madx_reference_bunch_num_particles=madx_reference_bunch_num_particles)
# generate_mad_bb_info(bb_df_b4, mode='from_dataframe',
#         madx_reference_bunch_num_particles=madx_reference_bunch_num_particles)

bb_dfs = {
    'b1': bb_df_b1,
    'b2': bb_df_b2,
    'b3': bb_df_b3,
    'b4': bb_df_b4}

if abs(z_crab_twiss)>0:
    crab_kicker_dict = crabbing_strong_beam(mad, bb_dfs,
            z_crab_twiss=z_crab_twiss,
            save_crab_twiss=True)
else:
    print('Crabbing of strong beam skipped!')

#if remove_dummy_lenses:
#    for beam in ['b1', 'b2']:
#        bbdf = bb_dfs[beam]
#        mad.input(f'seqedit, sequence={"lhc"+beam};')
#        mad.input('flatten;')
#        for nn in bbdf.elementName.values:
#            print(f'remove, element={nn}')
#            mad.input(f'remove, element={nn}')
#        mad.input('flatten;')
#        mad.input(f'endedit;')

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
from pymask.beambeam import crabbing_strong_beam_xsuite
from pymask.beambeam import setup_beam_beam_in_line

with open('../hl_lhc_collisions_000_b1_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b1 = json.load(fid)
with open('../hl_lhc_collisions_001_b4_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b4 = json.load(fid)
line_b1 = xt.Line.from_dict(dct_b1)
line_b4 = xt.Line.from_dict(dct_b4)


line_b1.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)
line_b4.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)


ip_names=['ip1', 'ip2', 'ip5', 'ip8']
num_long_range_elems_per_side=[25, 20, 25, 20]
harmonic_number=35640
bunch_spacing_buckets=10
num_slices_head_on=11
bunch_num_particles=2.2e11
sigmaz_m = 0.076
nemitt_x = 2.5e-6
nemitt_y = 2.5e-6
crab_strong_beam = True

circumference = line_b1.get_length()

from temp_module import install_beambeam_elements_in_lines

bb_df_b1_ret, bb_df_b2_ret = install_beambeam_elements_in_lines(line_b1, line_b4, ip_names,
            circumference, harmonic_number, bunch_spacing_buckets,
            num_long_range_elems_per_side, num_slices_head_on,
            bunch_num_particles, sigmaz_m)

tracker_b1 = line_b1.build_tracker()
tracker_b4 = line_b4.build_tracker()

keep_columns = ['beam', 'other_beam', 'ip_name', 'elementName', 'other_elementName', 'label',
                'self_num_particles', 'self_particle_charge', 'self_relativistic_beta',
                'identifier', 's_crab']
bb_df_b1 = bb_df_b1_ret[keep_columns].copy()
bb_df_b2 = bb_df_b2_ret[keep_columns].copy()

twiss_b1 = tracker_b1.twiss()
twiss_b4 = tracker_b4.twiss()
twiss_b2 = twiss_b4.reverse()

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
import pdb; pdb.set_trace()
for bb_df in [bb_df_b1, bb_df_b2]:
    compute_separations(bb_df)
    compute_dpx_dpy(bb_df)
    compute_local_crossing_angle_and_plane(bb_df)
    compute_xma_yma(bb_df)

# Get bb dataframe and mad model (with dummy bb) for beam 3 and 4
bb_df_b3 = get_counter_rotating(bb_df_b1)
bb_df_b4 = get_counter_rotating(bb_df_b2)

bb_dfs = {
    'b1': bb_df_b1,
    'b2': bb_df_b2,
    'b3': bb_df_b3,
    'b4': bb_df_b4}

if crab_strong_beam:
    crabbing_strong_beam_xsuite(bb_dfs,
        tracker_b1, tracker_b4)
else:
    print('Crabbing of strong beam skipped!')

setup_beam_beam_in_line(line_b1, bb_df_b1, bb_coupling=False)
setup_beam_beam_in_line(line_b4, bb_df_b4, bb_coupling=False)

xf.configure_orbit_dependent_parameters_for_bb(tracker=tracker_b1,
                    particle_on_co=twiss_b1.particle_on_co)
xf.configure_orbit_dependent_parameters_for_bb(tracker=tracker_b4,
                    particle_on_co=twiss_b4.particle_on_co)


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

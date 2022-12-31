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

bb_df_b1_ret, bb_df_b2_ret = install_beambeam_elements_in_lines(
            line_b1, line_b4, ip_names,
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


import json
import numpy as np
import xtrack as xt
import xpart as xp
from cpymad.madx import Madx

with open('../hl_lhc_collisions_000_b1_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b1 = json.load(fid)
with open('../hl_lhc_collisions_001_b4_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b4 = json.load(fid)
line_b1 = xt.Line.from_dict(dct_b1)
line_b4 = xt.Line.from_dict(dct_b4)

# mad_b1 = Madx()
# mad_b1.call('../hl_lhc_collisions_000_b1_no_bb/final_seq.madx')
# mad_b1.use('lhcb1')
# mad_b1.twiss()
# mad_b4 = Madx()
# mad_b4.call('../hl_lhc_collisions_001_b4_no_bb/final_seq.madx')
# mad_b4.use('lhcb2')
# mad_b4.twiss()
#line_b1 = xt.Line.from_madx_sequence(mad_b1.sequence.lhcb1, deferred_expressions=True)
#line_b4 = xt.Line.from_madx_sequence(mad_b4.sequence.lhcb2, deferred_expressions=True)


line_b1.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)
line_b4.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, p0c=7e12)

tracker_b1 = line_b1.build_tracker()
tracker_b4 = line_b4.build_tracker()

tw_b1 = tracker_b1.twiss()
tw_b4 = tracker_b4.twiss()

import pandas as pd
dfb2 = pd.read_parquet('../hl_lhc_collisions_001_b4_no_bb/twiss_b2_for_b4check_seq_lhcb2.parquet')

tw_b2 = tw_b4.mirror()

Ws = np.array(tw_b2.W_matrix)
alfx_from_wmat_b2 = -Ws[:, 0, 0] * Ws[:, 1, 0] - Ws[:, 0, 1] * Ws[:, 1, 1]

# Go to flat machine in xtrack
tracker_b1.vars['on_x1'] = 0.0
tracker_b1.vars['on_x2'] = 0.0
tracker_b1.vars['on_x5'] = 0.0
tracker_b1.vars['on_x8'] = 0.0
tracker_b1.vars['on_sep1'] = 0.0
tracker_b1.vars['on_sep2'] = 0.0
tracker_b1.vars['on_sep5'] = 0.0
tracker_b1.vars['on_sep8'] = 0.0
tracker_b1.vars['on_disp'] = 0.0
tracker_b1.vars['on_alice'] = 0
tracker_b1.vars['on_lhcb'] = 0 #  This one seems not to work
tracker_b1.vars['on_crab1'] = 0
tracker_b1.vars['on_crab5'] = 0
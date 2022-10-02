import json
import xtrack as xt
import xpart as xp

with open('../hl_lhc_collisions_000_b1_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b1 = json.load(fid)

with open('../hl_lhc_collisions_001_b4_no_bb/xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dct_b4 = json.load(fid)

line_b1 = xt.Line.from_dict(dct_b1)
line_b4 = xt.Line.from_dict(dct_b4)

line_b1.particle_ref = xp.Particles.from_dict(dct_b1['particle_on_tracker_co'])
line_b4.particle_ref = xp.Particles.from_dict(dct_b4['particle_on_tracker_co'])

tracker_b1 = line_b1.build_tracker()
tracker_b4 = line_b4.build_tracker()

tw_b1 = tracker_b1.twiss()
tw_b4 = tracker_b4.twiss()

import pandas as pd
dfb2 = pd.read_parquet('../hl_lhc_collisions_001_b4_no_bb/twiss_b2_for_b4check_seq_lhcb2.parquet')

tw_b2 = tw_b4.mirror()


# tracker_b1.vars['on_x1'] = 0.0
# tracker_b1.vars['on_x2'] = 0.0
# tracker_b1.vars['on_x5'] = 0.0
# tracker_b1.vars['on_x8'] = 0.0
# tracker_b1.vars['on_sep1'] = 0.0
# tracker_b1.vars['on_sep2'] = 0.0
# tracker_b1.vars['on_sep5'] = 0.0
# tracker_b1.vars['on_sep8'] = 0.0
# tracker_b1.vars['on_disp'] = 0.0
# tracker_b1.vars['on_alice'] = 0
#tracker_b1.vars['on_lhcb'] = 0 #  This one seems not to work

#twflat = tracker_b1.twiss()



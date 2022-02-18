import json
import numpy as np
import NAFFlib
import helpers as hp
import matplotlib.pyplot as plt

import xtrack as xt
import xpart as xp
import sixtracktools


displace_x = [-5e-4, 5e-4]
displace_y = [3e-4, -3e-4]
num_turns = 5


with open('../xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    dict_line_xtrack = json.load(fid)

line = xt.Line.from_dict(dict_line_xtrack)

partCO = xp.Particles.from_dict(dict_line_xtrack['particle_on_tracker_co'])

(x_tbt_sixtrack, px_tbt_sixtrack, y_tbt_sixtrack, py_tbt_sixtrack,
 zeta_tbt_sixtrack, delta_tbt_sixtrack, extra) = hp.track_particle_sixtrack(
        partCO=partCO, Dx_wrt_CO_m=np.array(displace_x),
        Dpx_wrt_CO_rad=0,
        Dy_wrt_CO_m=np.array(displace_y), Dpy_wrt_CO_rad=0.,
        Dzeta_wrt_CO_m=0., Ddelta_wrt_CO=0., n_turns=num_turns,
        input_folder='../')

tracker = xt.Tracker(line=line)

particles = xp.build_particles(particle_on_co=partCO,
        mode='shift',
        x=np.array(displace_x),
        y=np.array(displace_y))

tracker.track(particles, turn_by_turn_monitor=True, num_turns=num_turns)

print('Xtrack')
print(tracker.record_last_track.x.transpose())
print('Sixtrack')
print(x_tbt_sixtrack)

assert np.allclose(tracker.record_last_track.x[0, :], x_tbt_sixtrack[:,0],
       rtol=1e-15, atol=9e-11)
assert np.allclose(tracker.record_last_track.y[0, :], y_tbt_sixtrack[:,0],
       rtol=1e-15, atol=9e-11)
assert np.allclose(tracker.record_last_track.delta[0, :], delta_tbt_sixtrack[:,0],
       rtol=1e-15, atol=5e-11)

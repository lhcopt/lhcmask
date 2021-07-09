import json
import numpy as np
import NAFFlib
import helpers as hp
import footprint
import matplotlib.pyplot as plt

import xline as xl
import xtrack as xt
import sixtracktools


displace_x = [-5e-4, 5e-4]
displace_y = [3e-4, -3e-4]
num_turns = 5


with open('../xlines/line_bb_for_tracking.json', 'r') as fid:
    dict_line_xtrack = json.load(fid)

#with open('../xline/line_bb_dipole_cancelled.json', 'r') as fid:
#    dict_line_old = json.load(fid)

line = xl.Line.from_dict(dict_line_xtrack)
#line = xl.Line.from_dict(dict_line_old)

partCO = xl.Particles.from_dict(dict_line_xtrack['particle_on_tracker_co'])

(x_tbt_sixtrack, px_tbt_sixtrack, y_tbt_sixtrack, py_tbt_sixtrack,
 sigma_tbt_sixtrack, delta_tbt_sixtrack, extra) = hp.track_particle_sixtrack(
        partCO=partCO, Dx_wrt_CO_m=np.array(displace_x),
        Dpx_wrt_CO_rad=0,
        Dy_wrt_CO_m=np.array(displace_y), Dpy_wrt_CO_rad=0.,
        Dsigma_wrt_CO_m=0., Ddelta_wrt_CO=0., n_turns=num_turns,
        input_folder='../')

tracker = xt.Tracker(sequence=line)

part_track = partCO.copy()
part_track.x += np.array(displace_x)
part_track.y += np.array(displace_y)

particles = xt.Particles(**part_track.to_dict())
particles.particle_id = np.arange(particles.num_particles)

tracker.track(particles, turn_by_turn_monitor=True, num_turns=num_turns)

assert np.allclose(tracker.record_last_track.x[0, :], x_tbt_sixtrack[:,0],
       rtol=1e-15, atol=9e-11)
assert np.allclose(tracker.record_last_track.y[0, :], y_tbt_sixtrack[:,0],
       rtol=1e-15, atol=9e-11)
assert np.allclose(tracker.record_last_track.delta[0, :], delta_tbt_sixtrack[:,0],
       rtol=1e-15, atol=5e-11)

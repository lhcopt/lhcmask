import json
import numpy as np

import xtrack as xt
import xpart as xp

path = '../xsuite_lines/line_bb_for_tracking.json'

with open(path, 'r') as fid:
    line = xt.Line.from_dict(json.load(fid))

line._var_management = None # Remove knobs

line.remove_inactive_multipoles(inplace=True)
line.remove_zero_length_drifts(inplace=True)
line.merge_consecutive_drifts(inplace=True)
line.merge_consecutive_multipoles(inplace=True)

tracker = xt.Tracker(line=line,
                     extra_headers=['#define XTRACK_MULTIPOLE_NO_SYNRAD'],
                     )

num_turns = 10
num_particles = 860

particles = xp.Particles(x=np.linspace(-1e-4, 1e-4, num_particles), p0c=7e12)

import time
t1 = time.time()
tracker.track(particles, num_turns=num_turns)
t2 = time.time()

print(f'Elapsed time: {t2-t1:.2e} s; {(t2-t1)/num_particles/num_turns*1e6} us/part/turn')

print(f'{particles.x[3]=}')
assert np.all(particles.state > 0)

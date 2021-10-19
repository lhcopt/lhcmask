import json

import numpy as np

import xline as xl
import xtrack as xt

with open('../xlines/line_bb_for_tracking.json', 'r') as fid:
    line_dict = json.load(fid)

line = xl.Line.from_dict(line_dict)

partCO = xl.Particles.from_dict(line_dict['particle_on_tracker_co'])

tracker = xt.Tracker(sequence=line)
particles = xt.Particles(**partCO.to_dict())

for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker.track(particles)

WW = np.array(line_dict['WW_finite_diffs'])
WWinv = np.array(line_dict['WWInv_finite_diffs'])

assert np.max(np.abs(np.dot(WW, WWinv) - np.eye(6)))<1e-10


ampl_sigmas = 1.
norm_emit_x = 2.5e-6
geom_emit_x = norm_emit_x / particles.beta0 / particles.gamma0

n_part = 100
theta = np.linspace(0, 2*np.pi, n_part)
x_norm = ampl_sigmas * np.sqrt(geom_emit_x) * np.cos(theta)
px_norm = ampl_sigmas * np.sqrt(geom_emit_x) * np.sin(theta)

XX_norm= np.array([x_norm,
                   px_norm,
                   np.zeros(n_part),
                   np.zeros(n_part),
                   np.zeros(n_part),
                   np.zeros(n_part)])

XX = np.dot(WW, XX_norm)
pp = partCO.copy()

pp.x += XX[0, :]
pp.px += XX[1, :]
pp.y += XX[2, :]
pp.py += XX[3, :]
pp.zeta += XX[4, :]
pp.delta += XX[5, :]
particles_matched = xt.Particles(**pp.to_dict())


import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(particles_matched.x, particles_matched.px)

tracker.track(particles_matched, num_turns=10)

plt.plot(particles_matched.x, particles_matched.px)

plt.show()

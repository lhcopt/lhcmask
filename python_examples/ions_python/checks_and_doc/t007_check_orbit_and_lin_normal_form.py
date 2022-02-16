import json

import numpy as np

import xtrack as xt
import xpart as xp

with open('../xsuite_lines/line_bb_for_tracking.json', 'r') as fid:
    line_dict = json.load(fid)

line = xt.Line.from_dict(line_dict)

partCO = xp.Particles.from_dict(line_dict['particle_on_tracker_co'])

tracker = xt.Tracker(line=line)
particles = xp.Particles(**partCO.to_dict())

for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker.track(particles)

for nn in 'x px y py zeta delta'.split():
    assert np.abs(getattr(particles, nn) - getattr(partCO, nn)) < 3e-10

WW = np.array(line_dict['WW_finite_diffs'])
WWinv = np.array(line_dict['WWInv_finite_diffs'])
assert np.max(np.abs(np.dot(WW, WWinv) - np.eye(6)))<1e-10


ampl_sigmas = 0.2
norm_emit_x = 2.5e-6
geom_emit_x = norm_emit_x / particles.beta0 / particles.gamma0

n_part = 100
theta = np.linspace(0, 2*np.pi, n_part)
x_norm = ampl_sigmas * np.sqrt(geom_emit_x) * np.cos(theta)
px_norm = ampl_sigmas * np.sqrt(geom_emit_x) * np.sin(theta)

particles_matched  = xp.build_particles(particle_on_co=partCO,
                                        x_norm=x_norm, px_norm=px_norm,
                                        R_matrix=np.array(line_dict['RR_finite_diffs']))
particles_test = particles_matched.copy()
tracker.track(particles_test, num_turns=10)

i_matched = np.argmax(particles_matched.x)
i_test = np.argmax(particles_test.x)
assert np.abs(particles_test.x[i_test] - particles_matched.x[i_matched]) < 1e-6
assert np.abs(particles_test.px[i_test] - particles_matched.px[i_matched]) < 1e-7


import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(particles_matched.x, particles_matched.px)
plt.plot(particles_test.x, particles_test.px)

plt.show()

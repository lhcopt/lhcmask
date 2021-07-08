import sys
import numpy as np
import xline as xl
import xtrack as xt
from scipy.optimize import fsolve
sys.path.append('./modules')
import pymask as pm
Madx = pm.Madxp
mad = Madx(command_log="mad_final.log")
mad.call("final_seq.madx")
mad.use(sequence="lhcb1")
mad.twiss()
mad.readtable(file="final_errors.tfs", table="errtab")
mad.seterr(table="errtab")
mad.set(format=".15g")
mad.twiss(rmatrix = True)

separation_given_wrt_closed_orbit_4D = True
tw = mad.table['twiss']

seq_name = 'lhcb1'
mad_beam =  mad.sequence[seq_name].beam

assert mad_beam.deltap == 0, "Not implemented."

line = xl.Line.from_json(
        'xline/line_bb_dipole_not_cancelled.json')

particle_on_madx_co = xl.Particles(
    p0c = mad_beam.pc*1e9,
    q0 = mad_beam.charge,
    mass0 = mad_beam.mass*1e9,
    s = 0,
    x = tw.x[0],
    px = tw.px[0],
    y = tw.y[0],
    py = tw.py[0],
    tau = tw.t[0],
    ptau = tw.pt[0],
)

RR_madx = np.zeros([6,6])

for ii in range(6):
    for jj in range(6):
        RR_madx[ii, jj] = getattr(tw, f're{ii+1}{jj+1}')[0]


tracker = xt.Tracker(sequence=line)

for ee in tracker.line.elements:
    if ee.__class__.__name__.startswith('BeamBeam'):
         ee._temp_q0 = ee.q0
         ee.q0 = 0

def one_turn_map(p, particle_on_madx_co, tracker):
    xl_part = particle_on_madx_co.copy()
    xl_part.x = p[0]
    xl_part.px = p[1]
    xl_part.y = p[2]
    xl_part.py = p[3]
    xl_part.zeta = p[4]
    xl_part.delta = p[5]

    part = xt.Particles(**xl_part.to_dict())
    tracker.track(part)
    p_res = np.array([
           part.x[0],
           part.px[0],
           part.y[0],
           part.py[0],
           part.zeta[0],
           part.delta[0]])
    return p_res


print('Start CO search')
res = fsolve(lambda p: p - one_turn_map(p, particle_on_madx_co, tracker), x0=np.array(6*[0]))
print('Done CO search')

particle_on_co = particle_on_madx_co.copy()
particle_on_co.x = res[0]
particle_on_co.px = res[1]
particle_on_co.y = res[2]
particle_on_co.py = res[3]
particle_on_co.zeta = res[4]
particle_on_co.delta = res[5]

temp_particles = xt.Particles(**particle_on_co.to_dict())
for ii, ee in enumerate(tracker.line.elements):
    if ee.__class__.__name__ == 'BeamBeamBiGaussian2D':
          px_0 = temp_particles.px[0]
          py_0 = temp_particles.py[0]
          ee.q0 = ee._temp_q0

          if separation_given_wrt_closed_orbit_4D:
              ee.mean_x += temp_particles.x[0]
              ee.mean_y += temp_particles.y[0]
              line.elements[ii].x_bb = ee.mean_x
              line.elements[ii].y_bb = ee.mean_y

          ee.track(temp_particles)

          ee.d_px = temp_particles.px - px_0
          ee.d_py = temp_particles.py - py_0
          line.elements[ii].d_px = ee.d_px
          line.elements[ii].d_py = ee.d_py

          temp_particles.px -= ee.d_px
          temp_particles.py -= ee.d_py

    elif ee.__class__.__name__ == 'BeamBeamBiGaussian3D':
        ee.q0 = ee._temp_q0
        ee.x_CO = temp_particles.x[0]
        ee.px_CO = temp_particles.px[0]
        ee.y_CO = temp_particles.y[0]
        ee.py_CO = temp_particles.py[0]
        ee.sigma_CO = temp_particles.zeta[0]
        ee.delta_CO = temp_particles.delta[0]

        ee.track(temp_particles)

        ee.Dx_sub = temp_particles.x[0] - ee.x_CO
        ee.Dpx_sub = temp_particles.px[0] - ee.px_CO
        ee.Dy_sub = temp_particles.y[0] - ee.y_CO
        ee.Dpy_sub = temp_particles.py[0] - ee.py_CO
        ee.Dsigma_sub = temp_particles.zeta[0] - ee.sigma_CO
        ee.Ddelta_sub = temp_particles.delta[0] - ee.delta_CO

        temp_particles.x[0] = ee.x_CO
        temp_particles.px[0] = ee.px_CO
        temp_particles.y[0] = ee.y_CO
        temp_particles.py[0] = ee.py_CO
        temp_particles.zeta[0] = ee.sigma_CO
        temp_particles.delta[0] = ee.delta_CO

        line.elements[ii].x_co = ee.x_CO
        line.elements[ii].px_co = ee.px_CO
        line.elements[ii].y_co = ee.y_CO
        line.elements[ii].py_co = ee.py_CO
        line.elements[ii].zeta_co = ee.sigma_CO
        line.elements[ii].delta_co = ee.delta_CO

        line.elements[ii].d_x = ee.Dx_sub
        line.elements[ii].d_px = ee.Dpx_sub
        line.elements[ii].d_y = ee.Dy_sub
        line.elements[ii].d_py = ee.Dpy_sub
        line.elements[ii].d_zeta = ee.Dsigma_sub
        line.elements[ii].d_delta = ee.Ddelta_sub
    else:
        ee.track(temp_particles)


dict_line_xtrack = line.to_dict()

dict_line_xtrack['particle_on_co'] = particle_on_co.to_dict()

import pickle
with open('line_xtrack.pkl', 'wb') as fid:
    pickle.dump(dict_line_xtrack, fid)











# checks

particles = xt.Particles(**particle_on_co.to_dict())

print('Test closed orbit')
for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker.track(particles)

tracker_line_xtrack = xt.Tracker(sequence=xl.Line.from_dict(dict_line_xtrack)) 
particles = xt.Particles(**particle_on_co.to_dict())
print('Test closed orbit')
for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker_line_xtrack.track(particles)

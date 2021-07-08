import json

import numpy as np

import xline as xl
import xtrack as xt

line = xl.Line.from_json('../xline/line_bb_dipole_cancelled.json')

with open('../optics_orbit_at_start_ring.json', 'r') as fid:
    ddd = json.load(fid)
ddd['p0c'] =  ddd['p0c_eV']
partCO = xl.Particles.from_dict(ddd)

tracker = xt.Tracker(sequence=line)
particles = xt.Particles(**partCO.to_dict())

for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker.track(particles)


p0 = np.array([
           particles.x[0],
           particles.px[0],
           particles.y[0],
           particles.py[0],
           particles.zeta[0],
           particles.delta[0]])
def f(p):
    part = xt.Particles(
            p0c = partCO.p0c,
            x = p[0],
            px = p[1],
            y = p[2],
            py = p[3],
            zeta = p[4],
            delta = p[5])
    tracker.track(part)
    p_res = np.array([
           part.x[0],
           part.px[0],
           part.y[0],
           part.py[0],
           part.zeta[0],
           part.delta[0]])
    return p - p_res

from scipy.optimize import fsolve

res = fsolve(f, p0)
particles = xt.Particles(
            p0c = partCO.p0c,
            x = res[0],
            px = res[1],
            y = res[2],
            py = res[3],
            zeta = res[4],
            delta = res[5])
for _ in range(10):
    print(particles.at_turn[0], particles.x[0], particles.y[0],
          particles.zeta[0])
    tracker.track(particles)

p0 = res
II = np.eye(6)
JJ = np.array(
        [[ 0,-1, 0, 0, 0, 0],
         [ 1, 0, 0, 0, 0, 0],
         [ 0, 0, 0,-1, 0, 0],
         [ 0, 0, 1, 0, 0, 0],
         [ 0, 0, 0, 0, 0, -1],
         [ 0, 0, 0, 0, 1, 0]])

def f2(p):
    part = xt.Particles(
            p0c = partCO.p0c,
            x = p[0],
            px = p[1],
            y = p[2],
            py = p[3],
            zeta = p[4],
            delta = p[5])
    tracker.track(part)
    p_res = np.array([
           part.x[0],
           part.px[0],
           part.y[0],
           part.py[0],
           part.zeta[0],
           part.delta[0]])
    return p_res
RR = np.zeros((6, 6), dtype=np.float64)
for jj,dd in enumerate([1e-9,1e-12,1e-9,1e-12,1e-9,1e-9]):
    pd=p0+II[jj]*dd
    RR[:,jj]=(f2(pd)-p0)/dd


import json
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

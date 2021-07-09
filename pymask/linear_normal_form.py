import numpy as np
from scipy.optimize import fsolve

def _one_turn_map(p, particle_on_madx_co, tracker):
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

def _find_closed_orbit_from_tracker(tracker, particle_co_guess_dict):
    particle_co_guess = xl.Particles.from_dict(particle_co_guess_dict)
    print('Start CO search')
    res = fsolve(lambda p: p - _one_turn_map(p, particle_co_guess, tracker),
          x0=np.array([particle_co_guess.x, particle_co_guess.px,
                       particle_co_guess.y, particle_co_guess.py,
                       particle_co_guess.zeta, particle_co_guess.delta]))
    print('Done CO search')
    particle_on_co = particle_co_guess.copy()
    particle_on_co.x = res[0]
    particle_on_co.px = res[1]
    particle_on_co.y = res[2]
    particle_on_co.py = res[3]
    particle_on_co.zeta = res[4]
    particle_on_co.delta = res[5]

    return particle_on_co

def compute_R_matrix_finite_differences(one_turn_map,
        particle_on_co, tracker,
        dx=1e-9, dpx=1e-12, dy=1e-9, dpy=1e-12,
        dzeta=1e-9, ddelta=1e-9):

    # Find R matrix
    p_0 = np.array([
           particle_on_co.x[0],
           particle_on_co.px[0],
           particle_on_co.y[0],
           particle_on_co.py[0],
           particle_on_co.zeta[0],
           particle_on_co.delta[0]])
    II = np.eye(6)
    RR = np.zeros((6, 6), dtype=np.float64)
    for jj, dd in enumerate([dx, dpx, dy, dpy, dzeta, ddelta]):
        pd=p0+II[jj]*dd
        RR[:,jj]=(one_turn_map(pd, particle_on_co, tracker)-p0)/dd
    return RR

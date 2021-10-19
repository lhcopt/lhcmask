import numpy as np
from scipy.optimize import fsolve

import xline as xl
import xtrack as xt

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

def find_closed_orbit_from_tracker(tracker, particle_co_guess_dict):
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

def compute_R_matrix_finite_differences(
        particle_on_co, tracker,
        dx=1e-9, dpx=1e-12, dy=1e-9, dpy=1e-12,
        dzeta=1e-9, ddelta=1e-9,
        symplectify=True):

    # Find R matrix
    p0 = np.array([
           particle_on_co.x,
           particle_on_co.px,
           particle_on_co.y,
           particle_on_co.py,
           particle_on_co.zeta,
           particle_on_co.delta])
    II = np.eye(6)
    RR = np.zeros((6, 6), dtype=np.float64)
    for jj, dd in enumerate([dx, dpx, dy, dpy, dzeta, ddelta]):
        RR[:,jj]=(_one_turn_map(p0+II[jj]*dd, particle_on_co, tracker)-
                  _one_turn_map(p0-II[jj]*dd, particle_on_co, tracker))/(2*dd)

    if not 0.999 < np.linalg.det(RR) < 1.001:
        raise ValueError('The determinant of the RR is out tolerance.')
    
    if symplectify:
        return healy_symplectify(RR)
    else:
        return RR

def healy_symplectify(M):
    # https://accelconf.web.cern.ch/e06/PAPERS/WEPCH152.PDF
    print("Symplectifying linear One-Turn-Map...")

    print("Before symplectifying: det(M) = {}".format(np.linalg.det(M)))
    I = np.identity(6)

    S = np.array(
        [
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
        ]
    )

    V = np.matmul(S, np.matmul(I - M, np.linalg.inv(I + M)))
    W = (V + V.T) / 2
    if np.linalg.det(I - np.matmul(S, W)) != 0:
        M_new = np.matmul(I + np.matmul(S, W),
                          np.linalg.inv(I - np.matmul(S, W)))
    else:
        print("WARNING: det(I - SW) = 0!")
        V_else = np.matmul(S, np.matmul(I + M, np.linalg.inv(I - M)))
        W_else = (V_else + V_else.T) / 2
        M_new = -np.matmul(
            I + np.matmul(S, W_else), np.linalg(I - np.matmul(S, W_else))
        )

    print("After symplectifying: det(M) = {}".format(np.linalg.det(M_new)))
    return M_new

S = np.array([[0., 1., 0., 0., 0., 0.],
              [-1., 0., 0., 0., 0., 0.],
              [ 0., 0., 0., 1., 0., 0.],
              [ 0., 0.,-1., 0., 0., 0.],
              [ 0., 0., 0., 0., 0., 1.],
              [ 0., 0., 0., 0.,-1., 0.]])

######################################################
### Implement Normalization of fully coupled motion ##
######################################################

def Rot2D(mu):
    return np.array([[ np.cos(mu), np.sin(mu)],
                     [-np.sin(mu), np.cos(mu)]])

def compute_linear_normal_form(M):
    w0, v0 = np.linalg.eig(M)

    a0 = np.real(v0)
    b0 = np.imag(v0)

    index_list = [0,5,1,2,3,4]

    ##### Sort modes in pairs of conjugate modes #####

    conj_modes = np.zeros([3,2], dtype=np.int64)
    for j in [0,1]:
        conj_modes[j,0] = index_list[0]
        del index_list[0]

        min_index = 0
        min_diff = abs(np.imag(w0[conj_modes[j,0]] + w0[index_list[min_index]]))
        for i in range(1,len(index_list)):
            diff = abs(np.imag(w0[conj_modes[j,0]] + w0[index_list[i]]))
            if min_diff > diff:
                min_diff = diff
                min_index = i

        conj_modes[j,1] = index_list[min_index]
        del index_list[min_index]

    conj_modes[2,0] = index_list[0]
    conj_modes[2,1] = index_list[1]

    ##################################################
    #### Select mode from pairs with positive (real @ S @ imag) #####

    modes = np.empty(3, dtype=np.int64)
    for ii,ind in enumerate(conj_modes):
        if np.matmul(np.matmul(a0[:,ind[0]], S), b0[:,ind[0]]) > 0:
            modes[ii] = ind[0]
        else:
            modes[ii] = ind[1]

    ##################################################
    #### Sort modes such that (1,2,3) is close to (x,y,zeta) ####
    for i in [1,2]:
        if abs(v0[:,modes[0]])[0] < abs(v0[:,modes[i]])[0]:
            modes[0], modes[i] = modes[i], modes[0]

    if abs(v0[:,modes[1]])[2] < abs(v0[:,modes[2]])[2]:
        modes[2], modes[1] = modes[1], modes[2]

    ##################################################
    #### Rotate eigenvectors to the Courant-Snyder parameterization ####
    phase0 = np.log(v0[0,modes[0]]).imag
    phase1 = np.log(v0[2,modes[1]]).imag
    phase2 = np.log(v0[4,modes[2]]).imag

    v0[:,modes[0]] *= np.exp(-1.j*phase0)
    v0[:,modes[1]] *= np.exp(-1.j*phase1)
    v0[:,modes[2]] *= np.exp(-1.j*phase2)

    ##################################################
    #### Construct W #################################

    a1 = v0[:,modes[0]].real
    a2 = v0[:,modes[1]].real
    a3 = v0[:,modes[2]].real
    b1 = v0[:,modes[0]].imag
    b2 = v0[:,modes[1]].imag
    b3 = v0[:,modes[2]].imag

    n1 = 1./np.sqrt(np.matmul(np.matmul(a1, S), b1))
    n2 = 1./np.sqrt(np.matmul(np.matmul(a2, S), b2))
    n3 = 1./np.sqrt(np.matmul(np.matmul(a3, S), b3))

    a1 *= n1
    a2 *= n2
    a3 *= n3

    b1 *= n1
    b2 *= n2
    b3 *= n3

    W = np.array([a1,b1,a2,b2,a3,b3]).T
    W[abs(W) < 1.e-14] = 0. # Set very small numbers to zero.
    invW = np.matmul(np.matmul(S.T, W.T), S)

    ##################################################
    #### Get tunes and rotation matrix in the normalized coordinates ####

    mu1 = np.log(w0[modes[0]]).imag
    mu2 = np.log(w0[modes[1]]).imag
    mu3 = np.log(w0[modes[2]]).imag

    #q1 = mu1/(2.*np.pi)
    #q2 = mu2/(2.*np.pi)
    #q3 = mu3/(2.*np.pi)

    R = np.zeros_like(W)
    R[0:2,0:2] = Rot2D(mu1)
    R[2:4,2:4] = Rot2D(mu2)
    R[4:6,4:6] = Rot2D(mu3)
    ##################################################

    return W, invW, R

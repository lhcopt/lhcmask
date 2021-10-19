import numpy as np
import os
import sixtracktools
import xline
import time


def vectorize_all_coords(Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                         Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                         Dsigma_wrt_CO_m, Ddelta_wrt_CO):

    n_part = 1
    for var in [Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                Dsigma_wrt_CO_m, Ddelta_wrt_CO]:
        if hasattr(var, '__iter__'):
            if n_part == 1:
                n_part = len(var)
            assert len(var) == n_part

    Dx_wrt_CO_m = Dx_wrt_CO_m + np.zeros(n_part)
    Dpx_wrt_CO_rad = Dpx_wrt_CO_rad + np.zeros(n_part)
    Dy_wrt_CO_m = Dy_wrt_CO_m + np.zeros(n_part)
    Dpy_wrt_CO_rad = Dpy_wrt_CO_rad + np.zeros(n_part)
    Dsigma_wrt_CO_m = Dsigma_wrt_CO_m + np.zeros(n_part)
    Ddelta_wrt_CO = Ddelta_wrt_CO + np.zeros(n_part)

    return Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
     Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
     Dsigma_wrt_CO_m, Ddelta_wrt_CO


def track_particle_sixtrack(
                            partCO, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                            Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                            Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns,
                            input_folder='./'):

    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
    Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
    Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                         Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                         Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                         Dsigma_wrt_CO_m, Ddelta_wrt_CO)

    n_part = len(Dx_wrt_CO_m)

    wfold = 'temp_trackfun'

    if not os.path.exists(wfold):
        os.mkdir(wfold)

    os.system(f'cp {input_folder}/fort.* {wfold}')

    with open(f'{wfold}/fort.3', 'r') as fid:
        lines_f3 = fid.readlines()

    # Set initial coordinates
    i_start_ini = None
    for ii, ll in enumerate(lines_f3):
        if ll.startswith('INITIAL COO'):
            i_start_ini = ii
            break

    lines_f3[i_start_ini + 2] = '    0.\n'
    lines_f3[i_start_ini + 3] = '    0.\n'
    lines_f3[i_start_ini + 4] = '    0.\n'
    lines_f3[i_start_ini + 5] = '    0.\n'
    lines_f3[i_start_ini + 6] = '    0.\n'
    lines_f3[i_start_ini + 7] = '    0.\n'

    lines_f3[i_start_ini + 2 + 6] = '    0.\n'
    lines_f3[i_start_ini + 3 + 6] = '    0.\n'
    lines_f3[i_start_ini + 4 + 6] = '    0.\n'
    lines_f3[i_start_ini + 5 + 6] = '    0.\n'
    lines_f3[i_start_ini + 6 + 6] = '    0.\n'
    lines_f3[i_start_ini + 7 + 6] = '    0.\n'

    lines_f13 = []

    for i_part in range(n_part):
        if type(partCO) is xline.Particles:
            temp_part = partCO.copy()
        else:
            temp_part = xline.Particles(**partCO)
        temp_part.x     += Dx_wrt_CO_m[i_part]
        temp_part.px    += Dpx_wrt_CO_rad[i_part]
        temp_part.y     += Dy_wrt_CO_m[i_part]
        temp_part.py    += Dpy_wrt_CO_rad[i_part]
        temp_part.sigma += Dsigma_wrt_CO_m[i_part]
        temp_part.delta += Ddelta_wrt_CO[i_part]

        lines_f13.append('%.10e\n' % ((temp_part.x) * 1e3))
        lines_f13.append('%.10e\n' % ((temp_part.px) * temp_part.rpp * 1e3))
        lines_f13.append('%.10e\n' % ((temp_part.y) * 1e3))
        lines_f13.append('%.10e\n' % ((temp_part.py) * temp_part.rpp * 1e3))
        lines_f13.append('%.10e\n' % ((temp_part.sigma) * 1e3))
        lines_f13.append('%.10e\n' % ((temp_part.delta)))
        if i_part % 2 == 1:
            lines_f13.append('%.10e\n' % (temp_part.energy0 * 1e-6))
            lines_f13.append('%.10e\n' % (prev_part.energy * 1e-6))
            lines_f13.append('%.10e\n' % (temp_part.energy * 1e-6))
        prev_part = temp_part

    with open(wfold + '/fort.13', 'w') as fid:
        fid.writelines(lines_f13)

    if np.mod(n_part, 2) != 0:
        raise ValueError('SixTrack does not like this!')

    i_start_tk = None
    for ii, ll in enumerate(lines_f3):
        if ll.startswith('TRACKING PAR'):
            i_start_tk = ii
            break
    # Set number of turns and number of particles
    temp_list = lines_f3[i_start_tk + 1].split(' ')
    temp_list[0] = '%d' % n_turns
    temp_list[2] = '%d' % (n_part / 2)
    lines_f3[i_start_tk + 1] = ' '.join(temp_list)
    # Set number of idfor = 2
    temp_list = lines_f3[i_start_tk + 2].split(' ')
    temp_list[2] = '2'
    lines_f3[i_start_tk + 2] = ' '.join(temp_list)

    # Setup turn-by-turn dump
    i_start_dp = None
    for ii, ll in enumerate(lines_f3):
        if ll.startswith('DUMP'):
            i_start_dp = ii
            break

    lines_f3[i_start_dp + 1] = 'StartDUMP 1 664 101 dumtemp.dat\n'

    with open(wfold + '/fort.3', 'w') as fid:
        fid.writelines(lines_f3)

    os.system('(cd temp_trackfun; sixtrack >fort.6)')

    # Load sixtrack tracking data
    sixdump_all = sixtracktools.SixDump101('%s/dumtemp.dat' % wfold)

    x_tbt = np.zeros((n_turns, n_part))
    px_tbt = np.zeros((n_turns, n_part))
    y_tbt = np.zeros((n_turns, n_part))
    py_tbt = np.zeros((n_turns, n_part))
    sigma_tbt = np.zeros((n_turns, n_part))
    delta_tbt = np.zeros((n_turns, n_part))

    for i_part in range(n_part):
        sixdump_part = sixdump_all[i_part::n_part]
        x_tbt[:, i_part] = sixdump_part.x
        px_tbt[:, i_part] = sixdump_part.px
        y_tbt[:, i_part] = sixdump_part.y
        py_tbt[:, i_part] = sixdump_part.py
        sigma_tbt[:, i_part] = sixdump_part.sigma
        delta_tbt[:, i_part] = sixdump_part.delta

    extra = {}

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt, extra


def track_particle_xline(line, part, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                              Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                              Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns, verbose=False):

    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dsigma_wrt_CO_m, Ddelta_wrt_CO)

    part.x += Dx_wrt_CO_m
    part.px += Dpx_wrt_CO_rad
    part.y += Dy_wrt_CO_m
    part.py += Dpy_wrt_CO_rad
    part.sigma += Dsigma_wrt_CO_m
    part.delta += Ddelta_wrt_CO

    x_tbt = []
    px_tbt = []
    y_tbt = []
    py_tbt = []
    sigma_tbt = []
    delta_tbt = []

    for i_turn in range(n_turns):
        if verbose:
            print('Turn %d/%d' % (i_turn, n_turns))

        x_tbt.append(part.x.copy())
        px_tbt.append(part.px.copy())
        y_tbt.append(part.y.copy())
        py_tbt.append(part.py.copy())
        sigma_tbt.append(part.sigma.copy())
        delta_tbt.append(part.delta.copy())

        line.track(part)

    x_tbt = np.array(x_tbt)
    px_tbt = np.array(px_tbt)
    y_tbt = np.array(y_tbt)
    py_tbt = np.array(py_tbt)
    sigma_tbt = np.array(sigma_tbt)
    delta_tbt = np.array(delta_tbt)

    extra = {}

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt

def track_particle_xtrack(
                            line, partCO, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                            Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                            Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns,
                            _context=None):


    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dsigma_wrt_CO_m, Ddelta_wrt_CO)

    import xtrack as xt
    tracker = xt.Tracker(_context=_context, sequence=line)
    particles = xt.Particles(_context=_context,
            p0c = partCO.p0c,
            x = partCO.x + Dx_wrt_CO_m,
            px = partCO.px + Dpx_wrt_CO_rad,
            y = partCO.y + Dy_wrt_CO_m,
            py = partCO.py + Dpy_wrt_CO_rad,
            zeta = partCO.zeta + Dsigma_wrt_CO_m,
            delta = partCO.delta + Ddelta_wrt_CO)

    print('Start track')
    tracker.track(particles, num_turns=n_turns, turn_by_turn_monitor=True)
    print('Done track')

    #print(res.particles[0])
    x_tbt = tracker.record_last_track.x.copy().T
    px_tbt = tracker.record_last_track.px.copy().T
    y_tbt = tracker.record_last_track.y.copy().T
    py_tbt = tracker.record_last_track.py.copy().T
    sigma_tbt = tracker.record_last_track.zeta.copy().T
    delta_tbt = tracker.record_last_track.delta.copy().T

    print('Done loading!')

    extra = {'particles': particles, 'tracker': tracker}

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt, extra


def betafun_from_ellip(x_tbt, px_tbt):

    x_max = np.max(x_tbt)
    mask = np.logical_and(np.abs(x_tbt) < x_max / 5., px_tbt > 0)
    x_masked = x_tbt[mask]
    px_masked = px_tbt[mask]
    ind_sorted = np.argsort(x_masked)
    x_sorted = np.take(x_masked, ind_sorted)
    px_sorted = np.take(px_masked, ind_sorted)

    px_cut = np.interp(0, x_sorted, px_sorted)

    beta_x = x_max / px_cut

    return beta_x, x_max, px_cut

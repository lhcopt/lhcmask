import numpy as np
import os
import sixtracktools
import pysixtrack
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
        if type(partCO) is pysixtrack.Particles:
            temp_part = partCO.copy()
        else:
            temp_part = pysixtrack.Particles(**partCO)
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

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt


def track_particle_pysixtrack(line, part, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
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

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt

def track_particle_sixtracklib(
                            line, partCO, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                            Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                            Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns,
                            device=None):


    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dsigma_wrt_CO_m, Ddelta_wrt_CO)

    #if type(partCO) is pysixtrack.Particles:
    #    temp_part = partCO.copy()
    #else:
    #    temp_part = pysixtrack.Particles(**partCO)

    import sixtracklib
    elements=sixtracklib.Elements()
    elements.BeamMonitor(num_stores=n_turns)
    elements.append_line(line)

    n_part = len(Dx_wrt_CO_m)

    # Build PyST particle

    ps = sixtracklib.ParticlesSet()
    p = ps.Particles(num_particles=n_part)

    for i_part in range(n_part):

        if type(partCO) is pysixtrack.Particles:
            part = partCO.copy()
        else:
            part = pysixtrack.Particles(**partCO)
        part.x += Dx_wrt_CO_m[i_part]
        part.px += Dpx_wrt_CO_rad[i_part]
        part.y += Dy_wrt_CO_m[i_part]
        part.py += Dpy_wrt_CO_rad[i_part]
        part.sigma += Dsigma_wrt_CO_m[i_part]
        part.delta += Ddelta_wrt_CO[i_part]

        part.partid = i_part
        part.state = 1
        part.elemid = 0
        part.turn = 0

        p.from_pysixtrack(part, i_part)

    if device is None:
        job = sixtracklib.TrackJob(elements, ps)
    else:
        job = sixtracklib.TrackJob(elements, ps, device=device)

    job.track(n_turns)
    job.collect()

    res = job.output

    #print(res.particles[0])
    x_tbt = res.particles[0].x.reshape(n_turns, n_part)
    px_tbt = res.particles[0].px.reshape(n_turns, n_part)
    y_tbt = res.particles[0].y.reshape(n_turns, n_part)
    py_tbt = res.particles[0].py.reshape(n_turns, n_part)
    sigma_tbt = res.particles[0].sigma.reshape(n_turns, n_part)
    delta_tbt = res.particles[0].delta.reshape(n_turns, n_part)

    # For now data are saved at the end of the turn by STlib and at the beginning by the others
    #x_tbt[1:, :] = x_tbt[:-1, :]
    #px_tbt[1:, :] = px_tbt[:-1, :]
    #y_tbt[1:, :] = y_tbt[:-1, :]
    #py_tbt[1:, :] = py_tbt[:-1, :]
    #sigma_tbt[1:, :] = sigma_tbt[:-1, :]
    #delta_tbt[1:, :] = delta_tbt[:-1, :]
    #x_tbt[0, :] = p.x
    #px_tbt[0, :] = p.px
    #y_tbt[0, :] = p.y
    #py_tbt[0, :] = p.py
    #sigma_tbt[0, :] = p.sigma
    #delta_tbt[0, :] = p.delta

    print('Done loading!')

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt


def track_particle_sixtracklib_long(
                            line, partCO, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                            Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                            Dzeta_wrt_CO_m, Ddelta_wrt_CO, n_turns,
                            device=None):

    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dzeta_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dzeta_wrt_CO_m, Ddelta_wrt_CO)


    #if type(partCO) is pysixtrack.Particles:
    #    part = partCO.copy()
    #else:
    #    part = pysixtrack.Particles(**partCO)
    n_turns_tbt=1000
    n_turns_to_store=1000
    skip_turns=int(n_turns)//n_turns_to_store
    if skip_turns == 0:
        skip_turns = 1


    import sixtracklib
    elements=sixtracklib.Elements()
    #sixtracklib.append_beam_monitors_to_lattice(beam_elements_buffer=elements.cbuffer,
    #                                            until_turn_elem_by_elem=0,
    #                                            until_turn_turn_by_turn=n_turns_tbt,
    #                                            until_turn=n_turns,
    #                                            skip_turns=skip_turns
    #                                           )
    elements.BeamMonitor(num_stores=n_turns_tbt,start=0,skip=1,is_rolling=False)
    elements.BeamMonitor(num_stores=n_turns_to_store,start=0,skip=skip_turns,is_rolling=False)
    elements.BeamMonitor(num_stores=1,start=0,skip=1,is_rolling=True)
    print(elements.get_elements())
    #elements.BeamMonitor(num_stores=n_turns)
    #elements.BeamMonitor(num_stores=n_turns_to_store)
    elements.append_line(line)

    n_stores0=elements.get_elements()[0].num_stores
    n_stores1=elements.get_elements()[1].num_stores
    n_stores2=elements.get_elements()[2].num_stores
    n_part = len(Dx_wrt_CO_m)

    # Build PyST particle

    ps = sixtracklib.ParticlesSet()
    p = ps.Particles(num_particles=n_part)

    for i_part in range(n_part):

        if type(partCO) is pysixtrack.Particles:
            part = partCO.copy()
        else:
            part = pysixtrack.Particles(**partCO)
        part.x += Dx_wrt_CO_m[i_part]
        part.px += Dpx_wrt_CO_rad[i_part]
        part.y += Dy_wrt_CO_m[i_part]
        part.py += Dpy_wrt_CO_rad[i_part]
        part.zeta += Dzeta_wrt_CO_m[i_part]
        part.delta += Ddelta_wrt_CO[i_part]

        part.partid = i_part
        part.state = 1
        part.elemid = 0
        part.turn = 0

        p.from_pysixtrack(part, i_part)

    if device is None:
        job = sixtracklib.TrackJob(elements, ps)
    else:
        job = sixtracklib.TrackJob(elements, ps, device=device)

    start_tracking_time = time.time()
    job.track(n_turns)
    end_tracking_time = time.time()
    job.collect()
    end_collecting_time = time.time()
    res = job.output

    #print(res.particles[0])
    #print(res.particles[1])

    x_tbt_first       = res.particles[0].x.reshape(n_stores0,n_part)    
    px_tbt_first      = res.particles[0].px.reshape(n_stores0,n_part)    
    y_tbt_first       = res.particles[0].y.reshape(n_stores0,n_part)    
    py_tbt_first      = res.particles[0].py.reshape(n_stores0,n_part)    
    zeta_tbt_first    = res.particles[0].zeta.reshape(n_stores0,n_part)    
    delta_tbt_first   = res.particles[0].delta.reshape(n_stores0,n_part)    
    at_turn_tbt_first = res.particles[0].at_turn.reshape(n_stores0,n_part)    
    state_tbt_first   = res.particles[0].state.reshape(n_stores0,n_part)    

    x_skip       = res.particles[1].x.reshape(n_stores1,n_part)    
    px_skip      = res.particles[1].px.reshape(n_stores1,n_part)    
    y_skip       = res.particles[1].y.reshape(n_stores1,n_part)    
    py_skip      = res.particles[1].py.reshape(n_stores1,n_part)    
    zeta_skip    = res.particles[1].zeta.reshape(n_stores1,n_part)    
    delta_skip   = res.particles[1].delta.reshape(n_stores1,n_part)    
    at_turn_skip = res.particles[1].at_turn.reshape(n_stores1,n_part)    
    state_skip   = res.particles[1].state.reshape(n_stores1,n_part)    

    x_last       = res.particles[2].x.reshape(n_stores2,n_part)    
    px_last      = res.particles[2].px.reshape(n_stores2,n_part)    
    y_last       = res.particles[2].y.reshape(n_stores2,n_part)    
    py_last      = res.particles[2].py.reshape(n_stores2,n_part)    
    zeta_last    = res.particles[2].zeta.reshape(n_stores2,n_part)    
    delta_last   = res.particles[2].delta.reshape(n_stores2,n_part)    
    at_turn_last = res.particles[2].at_turn.reshape(n_stores2,n_part)    
    state_last   = res.particles[2].state.reshape(n_stores2,n_part)    

    output_dict = {'x_tbt_first' : x_tbt_first,
                   'px_tbt_first' : px_tbt_first,
                   'y_tbt_first' : y_tbt_first,
                   'py_tbt_first' : py_tbt_first,
                   'zeta_tbt_first' : zeta_tbt_first,
                   'delta_tbt_first' : delta_tbt_first,
                   'at_turn_tbt_first' : at_turn_tbt_first,
#                   'state_tbt_first' : state_tbt_first,
                   'x_skip'      : x_skip,      
                   'px_skip'     : px_skip,    
                   'y_skip'      : y_skip,      
                   'py_skip'     : py_skip,     
                   'zeta_skip'   : zeta_skip,   
                   'delta_skip'  : delta_skip,  
                   'at_turn_skip': at_turn_skip,
                   'state_skip'  : state_skip,  
                   'x_last' : x_last,
                   'px_last' : px_last,
                   'y_last' : y_last,
                   'py_last' : py_last,
                   'zeta_last' : zeta_last,
                   'delta_last' : delta_last,
                   'at_turn_last' : at_turn_last,
#                   'state_tbt_last' : state_tbt_last,
                   'tracking_time_mins' : (end_tracking_time - start_tracking_time)/60.,
                   'collecting_time_mins' : (end_collecting_time - end_tracking_time)/60.,
                  }


    print('Done loading!')
    return output_dict

def track_particle_sixtracklib_firstlast(
                            line, partCO, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                            Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                            Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns,
                            device=None):

    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dsigma_wrt_CO_m, Ddelta_wrt_CO)


    #if type(partCO) is pysixtrack.Particles:
    #    part = partCO.copy()
    #else:
    #    part = pysixtrack.Particles(**partCO)

    n_turns_to_store=1000
    n_turns_tbt=1000
    #skip_turns=1000


    import sixtracklib
    elements=sixtracklib.Elements()
    #sixtracklib.append_beam_monitors_to_lattice(beam_elements_buffer=elements.cbuffer,
    #                                            until_turn_elem_by_elem=0,
    #                                            until_turn_turn_by_turn=n_turns_tbt,
    #                                            until_turn=n_turns,
    #                                            skip_turns=skip_turns
    #                                           )
    elements.BeamMonitor(num_stores=n_turns_tbt,start=0,skip=1,is_rolling=False)
    elements.BeamMonitor(num_stores=n_turns_to_store,start=0,skip=1,is_rolling=True)
    print(elements.get_elements())
    #elements.BeamMonitor(num_stores=n_turns)
    #elements.BeamMonitor(num_stores=n_turns_to_store)
    elements.append_line(line)

    n_stores=elements.get_elements()[1].num_stores
    n_part = len(Dx_wrt_CO_m)

    # Build PyST particle

    ps = sixtracklib.ParticlesSet()
    p = ps.Particles(num_particles=n_part)

    for i_part in range(n_part):

        if type(partCO) is pysixtrack.Particles:
            part = partCO.copy()
        else:
            part = pysixtrack.Particles(**partCO)
        part.x += Dx_wrt_CO_m[i_part]
        part.px += Dpx_wrt_CO_rad[i_part]
        part.y += Dy_wrt_CO_m[i_part]
        part.py += Dpy_wrt_CO_rad[i_part]
        part.sigma += Dsigma_wrt_CO_m[i_part]
        part.delta += Ddelta_wrt_CO[i_part]

        part.partid = i_part
        part.state = 1
        part.elemid = 0
        part.turn = 0

        p.from_pysixtrack(part, i_part)

    if device is None:
        job = sixtracklib.TrackJob(elements, ps)
    else:
        job = sixtracklib.TrackJob(elements, ps, device=device)

    start_tracking_time = time.time()
    job.track(n_turns)
    end_tracking_time = time.time()
    job.collect()
    end_collecting_time = time.time()
    res = job.output

    print(res.particles[0])
    print(res.particles[1])

    x_tbt_first       = res.particles[0].x.reshape(n_turns_tbt,n_part)    
    px_tbt_first      = res.particles[0].px.reshape(n_turns_tbt,n_part)    
    y_tbt_first       = res.particles[0].y.reshape(n_turns_tbt,n_part)    
    py_tbt_first      = res.particles[0].py.reshape(n_turns_tbt,n_part)    
    zeta_tbt_first    = res.particles[0].zeta.reshape(n_turns_tbt,n_part)    
    delta_tbt_first   = res.particles[0].delta.reshape(n_turns_tbt,n_part)    
    at_turn_tbt_first = res.particles[0].at_turn.reshape(n_turns_tbt,n_part)    
    state_tbt_first   = res.particles[0].state.reshape(n_turns_tbt,n_part)    

    x_tbt_last       = res.particles[1].x.reshape(n_stores,n_part)    
    px_tbt_last      = res.particles[1].px.reshape(n_stores,n_part)    
    y_tbt_last       = res.particles[1].y.reshape(n_stores,n_part)    
    py_tbt_last      = res.particles[1].py.reshape(n_stores,n_part)    
    zeta_tbt_last    = res.particles[1].zeta.reshape(n_stores,n_part)    
    delta_tbt_last   = res.particles[1].delta.reshape(n_stores,n_part)    
    at_turn_tbt_last = res.particles[1].at_turn.reshape(n_stores,n_part)    
    state_tbt_last   = res.particles[1].state.reshape(n_stores,n_part)    

    output_dict = {'x_tbt_first' : x_tbt_first,
                   'px_tbt_first' : px_tbt_first,
                   'y_tbt_first' : y_tbt_first,
                   'py_tbt_first' : py_tbt_first,
                   'zeta_tbt_first' : zeta_tbt_first,
                   'delta_tbt_first' : delta_tbt_first,
                   'at_turn_tbt_first' : at_turn_tbt_first,
#                   'state_tbt_first' : state_tbt_first,
                   'x_tbt_last' : x_tbt_last,
                   'px_tbt_last' : px_tbt_last,
                   'y_tbt_last' : y_tbt_last,
                   'py_tbt_last' : py_tbt_last,
                   'zeta_tbt_last' : zeta_tbt_last,
                   'delta_tbt_last' : delta_tbt_last,
                   'at_turn_tbt_last' : at_turn_tbt_last,
#                   'state_tbt_last' : state_tbt_last,
                   'tracking_time_mins' : (end_tracking_time - start_tracking_time)/60.,
                   'collecting_time_mins' : (end_collecting_time - end_tracking_time)/60.,
                  }


    print('Done loading!')
    return output_dict

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

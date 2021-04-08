import os
import pickle

import numpy as np

from . import beambeam as bb

def make_links(links_dict, force=False):
    for kk in links_dict.keys():
        if force:
            if os.path.exists(kk):
                os.remove(kk)
        os.symlink(os.path.abspath(links_dict[kk]), kk)

def get_pymask_configuration(mode):

    if mode=='b1_without_bb':
        beam_to_configure = 1
        sequences_to_check = ['lhcb1', 'lhcb2']
        sequence_to_track = 'lhcb1'
        generate_b4_from_b2 = False
        track_from_b4_mad_instance = False
        enable_bb_python = False
        enable_bb_legacy = False
        force_disable_check_separations_at_ips = False
    elif mode=='b1_with_bb':
        beam_to_configure = 1
        sequences_to_check = ['lhcb1', 'lhcb2']
        sequence_to_track = 'lhcb1'
        generate_b4_from_b2 = False
        track_from_b4_mad_instance = False
        enable_bb_python = True
        enable_bb_legacy = False
        force_disable_check_separations_at_ips = False
    elif mode=='b1_with_bb_legacy_macros':
        beam_to_configure = 1
        sequences_to_check = ['lhcb1', 'lhcb2']
        sequence_to_track = 'lhcb1'
        generate_b4_from_b2 = False
        track_from_b4_mad_instance = False
        enable_bb_python = False
        enable_bb_legacy = True
        force_disable_check_separations_at_ips = False
    elif mode == 'b4_without_bb':
        beam_to_configure = 4
        sequences_to_check = ['lhcb2']
        sequence_to_track = 'lhcb2'
        generate_b4_from_b2 = False
        track_from_b4_mad_instance = False
        enable_bb_python = False
        enable_bb_legacy = False
        force_disable_check_separations_at_ips = True
    elif mode == 'b4_from_b2_without_bb':
        beam_to_configure = 1
        sequences_to_check = ['lhcb1', 'lhcb2']
        sequence_to_track = 'lhcb2'
        generate_b4_from_b2 = True
        track_from_b4_mad_instance = True
        enable_bb_python = False
        enable_bb_legacy = False
        force_disable_check_separations_at_ips = False
    elif mode == 'b4_from_b2_with_bb':
        beam_to_configure = 1
        sequences_to_check = ['lhcb1', 'lhcb2']
        sequence_to_track = 'lhcb2'
        generate_b4_from_b2 = True
        track_from_b4_mad_instance = True
        enable_bb_python = True
        enable_bb_legacy = False
        force_disable_check_separations_at_ips = False
    else:
        raise ValueError(f'Mode "{mode}" not recognized!')

    return (
        beam_to_configure,
        sequences_to_check,
        sequence_to_track,
        generate_b4_from_b2,
        track_from_b4_mad_instance,
        enable_bb_python,
        enable_bb_legacy,
        force_disable_check_separations_at_ips,
    )

def configure_b4_from_b2(mad_b4, mad_b2):
    var_dicts_b2 = mad_b2.get_variables_dicts()
    var_dicts_b4 = mad_b4.get_variables_dicts()

    b2_const=var_dicts_b2['constants']
    b4_const=var_dicts_b4['constants']
    for nn in b2_const.keys():
        if nn[0]=='_':
            print(f'The constant {nn} cannot be assigned!')
        else:
            if nn not in b4_const.keys():
                mad_b4.input(f'const {nn}={b2_const[nn]:.50e}')

    # %% INDEPENDENT
    b2_indep=var_dicts_b2['independent_variables']
    b4_indep=var_dicts_b4['independent_variables']
    for nn in b2_indep.keys():
        mad_b4.input(f'{nn}={b2_indep[nn]:.50e}')

    # %% DEPENDENT
    b2_dep=var_dicts_b2['dependent_variables_expr']
    b4_dep=var_dicts_b4['dependent_variables_expr']
    for nn in b2_dep.keys():
        mad_b4.input(f'{nn}:={str(b2_dep[nn])}')

    # bv_aux and my my lhcbeam need to be defined explicitly
    mad_b4.input(f'bv_aux=-1')
    mad_b4.input(f'mylhcbeam=4')

    # Attach beam
    mad_b4.use('lhcb2')
    beam_command = str(mad_b2.sequence['lhcb2'].beam)
    assert(', bv=-1.0' in beam_command)
    beam_command = beam_command.replace(', bv=-1.0', ', bv=1.0')
    mad_b4.input(beam_command)
    mad_b4.use('lhcb2')

    # %% CHECKS
    var_dicts_b2 = mad_b2.get_variables_dicts()
    var_dicts_b4 = mad_b4.get_variables_dicts()

    b2_const=var_dicts_b2['constants']
    b4_const=var_dicts_b4['constants']
    for nn in b4_const.keys():
        assert b2_const[nn] == b4_const[nn]

    for nn in b2_const.keys():
        if nn not in b4_const.keys():
            print(f'Warning: b2 const {nn}={b2_const[nn]} is not in b4.')

    b2_indep=var_dicts_b2['independent_variables']
    b4_indep=var_dicts_b4['independent_variables']
    for nn in b2_indep.keys():
        if str(nn) in 'bv_aux mylhcbeam':
            continue
        assert b4_indep[nn] == b2_indep[nn]

    for nn in b4_indep.keys():
        if nn not in b2_indep.keys():
            print(f'Warning: b4 indep {nn}={b4_indep[nn]} is not in b2.')

    b2_dep=var_dicts_b2['dependent_variables_expr']
    b4_dep=var_dicts_b4['dependent_variables_expr']
    for nn in b2_dep.keys():
        if str(nn) in 'bv_aux mylhcbeam':
            continue
        assert str(b4_dep[nn]) == str(b2_dep[nn])

    for nn in b4_dep.keys():
        if nn not in b2_dep.keys():
            print(f'Warning: b4 dep {nn}={str(b4_dep[nn])} is not in b2.')

def check_twiss_value(twiss_df, element_name, keyword, target, tol):
    assert abs(twiss_df.loc[element_name][keyword] - target) < tol,\
                f'Check not passed on {keyword} at {element_name}'

def check_twiss_against_madvars(checks, twiss_df, variable_dicts):
    for cc in checks:
        check_twiss_value(twiss_df,
            element_name=cc['element_name'],
            keyword=cc['keyword'],
            target=variable_dicts['all_variables_val'][cc['varname']],
            tol=cc['tol'])

def check_separation_value(twiss_df_b1, twiss_df_b2, element_name,
        plane, target, tol):
    assert plane in 'xy'
    val = (twiss_df_b2.loc[element_name, plane]
            - twiss_df_b1.loc[element_name, plane])
    assert abs(val - target) < tol,\
                f'Check not passed on {plane} separation at {element_name}'

def check_separations_against_madvars(checks, twiss_df_b1, twiss_df_b2, variables_dict):
    for cc in checks:
        tol = cc['tol']
        target = variables_dict['all_variables_val'][cc['varname']]*cc['scale_factor']
        check_separation_value(twiss_df_b1, twiss_df_b2, cc['element_name'],
                cc['plane'], target, cc['tol'])

def generate_sixtrack_input(mad, seq_name, bb_df, output_folder,
        reference_bunch_charge_sixtrack_ppb,
        emitnx_sixtrack_um,
        emitny_sixtrack_um,
        sigz_sixtrack_m,
        sige_sixtrack,
        ibeco_sixtrack,
        ibtyp_sixtrack,
        lhc_sixtrack,
        ibbc_sixtrack,
        radius_sixtrack_multip_conversion_mad,
        skip_mad_use=False):

    six_fol_name = output_folder
    os.makedirs(six_fol_name, exist_ok=True)

    os.system('rm fc.*')
    if not skip_mad_use:
        mad.use(seq_name)
    mad.twiss()
    mad.input(f'sixtrack, cavall, radius={radius_sixtrack_multip_conversion_mad}')
    os.system(f'mv fc.* {six_fol_name}')
    os.system(f'cp {six_fol_name}/fc.2 {six_fol_name}/fc.2.old')

    with open(six_fol_name + '/fc.2', 'r') as fid:
        fc2lines = fid.readlines()

    for ii, ll in enumerate(fc2lines):
        llfields = ll.split()
        try:
            if int(llfields[1]) == 20:
                newll = ' '.join([
                    llfields[0],
                    llfields[1]]
                    + (len(llfields)-2)* ['0.0']
                    +['\n'])
                fc2lines[ii] = newll
        except ValueError:
            pass # line does not have an integer in the second field
        except IndexError:
            pass # line has less than two fields

    with open(six_fol_name + '/fc.2', 'w') as fid:
        fid.writelines(fc2lines)

    # http://sixtrack.web.cern.ch/SixTrack/docs/user_full/manual.php#Ch6.S6

    if bb_df is not None:
        sxt_df_4d = bb_df[bb_df['label']=='bb_lr'].copy()
        if len(sxt_df_4d) > 0:
            sxt_df_4d['h-sep [mm]'] = -sxt_df_4d['separation_x']*1e3
            sxt_df_4d['v-sep [mm]'] = -sxt_df_4d['separation_y']*1e3
            sxt_df_4d['strength-ratio'] = (sxt_df_4d['other_charge_ppb']
                    / reference_bunch_charge_sixtrack_ppb) * mad.globals.on_bb_charge
            sxt_df_4d['4dSxx [mm*mm]'] = sxt_df_4d['other_Sigma_11']*1e6
            sxt_df_4d['4dSyy [mm*mm]'] = sxt_df_4d['other_Sigma_33']*1e6
            sxt_df_4d['4dSxy [mm*mm]'] = sxt_df_4d['other_Sigma_13']*1e6
            sxt_df_4d['fort3entry'] = sxt_df_4d.apply(lambda x: ' '.join([
                    f"{x.elementName}",
                    '0',
                    f"{x['4dSxx [mm*mm]']}",
                    f"{x['4dSyy [mm*mm]']}",
                    f"{x['h-sep [mm]']}",
                    f"{x['v-sep [mm]']}",
                    f"{x['strength-ratio']}",
                    # f"{x['4dSxy [mm*mm]']}" Not really used
                    ]), axis=1)


        sxt_df_6d = bb_df[bb_df['label']=='bb_ho'].copy()
        if len(sxt_df_6d) > 0:
            sxt_df_6d['h-sep [mm]'] = -sxt_df_6d['separation_x']*1e3
            sxt_df_6d['v-sep [mm]'] = -sxt_df_6d['separation_y']*1e3
            sxt_df_6d['phi [rad]'] = sxt_df_6d['phi']
            sxt_df_6d['alpha [rad]'] = sxt_df_6d['alpha']
            sxt_df_6d['strength-ratio'] = (sxt_df_6d['other_charge_ppb']/reference_bunch_charge_sixtrack_ppb)* mad.globals.on_bb_charge
            sxt_df_6d['Sxx [mm*mm]'] = sxt_df_6d['other_Sigma_11'] *1e6
            sxt_df_6d['Sxxp [mm*mrad]'] = sxt_df_6d['other_Sigma_12'] *1e6
            sxt_df_6d['Sxpxp [mrad*mrad]'] = sxt_df_6d['other_Sigma_22'] *1e6
            sxt_df_6d['Syy [mm*mm]'] = sxt_df_6d['other_Sigma_33'] *1e6
            sxt_df_6d['Syyp [mm*mrad]'] = sxt_df_6d['other_Sigma_34'] *1e6
            sxt_df_6d['Sypyp [mrad*mrad]'] = sxt_df_6d['other_Sigma_44'] *1e6
            sxt_df_6d['Sxy [mm*mm]'] = sxt_df_6d['other_Sigma_13'] *1e6
            sxt_df_6d['Sxyp [mm*mrad]'] = sxt_df_6d['other_Sigma_14'] *1e6
            sxt_df_6d['Sxpy [mrad*mm]'] = sxt_df_6d['other_Sigma_23'] *1e6
            sxt_df_6d['Sxpyp [mrad*mrad]'] = sxt_df_6d['other_Sigma_24'] *1e6
            sxt_df_6d['fort3entry'] = sxt_df_6d.apply(lambda x: ' '.join([
                    f"{x.elementName}",
                    '1',
                    f"{x['phi [rad]']}",
                    f"{x['alpha [rad]']}",
                    f"{x['h-sep [mm]']}",
                    f"{x['v-sep [mm]']}",
                    '\n'
                    f"{x['Sxx [mm*mm]']}",
                    f"{x['Sxxp [mm*mrad]']}",
                    f"{x['Sxpxp [mrad*mrad]']}",
                    f"{x['Syy [mm*mm]']}",
                    f"{x['Syyp [mm*mrad]']}",
                    '\n',
                    f"{x['Sypyp [mrad*mrad]']}",
                    f"{x['Sxy [mm*mm]']}",
                    f"{x['Sxyp [mm*mrad]']}",
                    f"{x['Sxpy [mrad*mm]']}",
                    f"{x['Sxpyp [mrad*mrad]']}",
                    f"{x['strength-ratio']}",
                    ]), axis=1)

        f3_common_settings = ' '.join([
            f"{reference_bunch_charge_sixtrack_ppb}",
            f"{emitnx_sixtrack_um}",
            f"{emitny_sixtrack_um}",
            f"{sigz_sixtrack_m}",
            f"{sige_sixtrack}",
            f"{ibeco_sixtrack}",
            f"{ibtyp_sixtrack}",
            f"{lhc_sixtrack}",
            f"{ibbc_sixtrack}",
            ])

        f3_string = '\n'.join([
            'BEAM',
            'EXPERT',
            f3_common_settings])

        f3_string += '\n'

        list_entries = []
        if len(sxt_df_6d) > 0:
            list_entries += list(sxt_df_6d['fort3entry'].values)
        if len(sxt_df_4d) > 0:
            list_entries += list(sxt_df_4d['fort3entry'].values)

        f3_string += '\n'.join(list_entries)

        f3_string += '\nNEXT\n'

        with open(six_fol_name + '/fc.3', 'a') as fid:
            fid.write(f3_string)


def get_optics_and_orbit_at_start_ring(mad, seq_name, with_bb_forces=False,
        skip_mad_use=False):

    initial_bb_state = None

    try:
        initial_bb_state = mad.globals.on_bb_switch
        mad.globals.on_bb_switch = {True: 1, False: 0}[with_bb_forces]
    except AttributeError:
        print('Warning! on_bb_switch not present')

    # Twiss and get closed-orbit
    if not skip_mad_use:
        mad.use(sequence=seq_name)
    twiss_table = mad.twiss()

    if initial_bb_state is not None:
        mad.globals.on_bb_switch = initial_bb_state

    beta0 = mad.sequence[seq_name].beam.beta
    gamma0 = mad.sequence[seq_name].beam.gamma
    p0c_eV = mad.sequence[seq_name].beam.pc*1.e9

    optics_at_start_ring = {
            'beta': beta0,
            'gamma' : gamma0,
            'p0c_eV': p0c_eV,
            'betx': twiss_table.betx[0],
            'bety': twiss_table.bety[0],
            'alfx': twiss_table.alfx[0],
            'alfy': twiss_table.alfy[0],
            'dx': twiss_table.dx[0],
            'dy': twiss_table.dy[0],
            'dpx': twiss_table.dpx[0],
            'dpy': twiss_table.dpy[0],
            'x' : twiss_table.x[0],
            'px' : twiss_table.px[0],
            'y' : twiss_table.y[0],
            'py' : twiss_table.py[0],
            't' : twiss_table.t[0],
            'pt' : twiss_table.pt[0],
            #convert tau, pt to sigma,delta
            'sigma' : beta0 * twiss_table.t[0],
            'delta' : ((twiss_table.pt[0]**2 +
                 2.*twiss_table.pt[0]/beta0) + 1.)**0.5 - 1.
            }
    return optics_at_start_ring

def generate_pysixtrack_line_with_bb(mad, seq_name, bb_df,
        closed_orbit_method='from_mad', pickle_lines_in_folder=None,
        skip_mad_use=False):

    opt_and_CO = get_optics_and_orbit_at_start_ring(mad, seq_name,
            with_bb_forces=False, skip_mad_use=True)

    # Build pysixtrack model
    import pysixtrack
    pysxt_line = pysixtrack.Line.from_madx_sequence(
        mad.sequence[seq_name])

    if bb_df is not None:
        bb.setup_beam_beam_in_line(pysxt_line, bb_df, bb_coupling=False)

    # Temporary fix due to bug in mad loader
    cavities, cav_names = pysxt_line.get_elements_of_type(
            pysixtrack.elements.Cavity)
    for cc, nn in zip(cavities, cav_names):
        if cc.frequency ==0.:
            ii_mad = mad.sequence[seq_name].element_names().index(nn)
            cc_mad = mad.sequence[seq_name].elements[ii_mad]
            f0_mad = mad.sequence[seq_name].beam.freq0 * 1e6 # mad has it in MHz
            cc.frequency = f0_mad*cc_mad.parent.harmon

    mad_CO = np.array([opt_and_CO[kk] for kk in ['x', 'px', 'y', 'py', 'sigma', 'delta']])

    pysxt_line.disable_beambeam()
    part_on_CO = pysxt_line.find_closed_orbit(
        guess=mad_CO, p0c=opt_and_CO['p0c_eV'],
        method={'from_mad': 'get_guess', 'from_tracking': 'Nelder-Mead'}[closed_orbit_method])
    pysxt_line.enable_beambeam()

    pysxt_line_bb_dipole_cancelled = pysxt_line.copy()

    pysxt_line_bb_dipole_cancelled.beambeam_store_closed_orbit_and_dipolar_kicks(
        part_on_CO,
        separation_given_wrt_closed_orbit_4D=True,
        separation_given_wrt_closed_orbit_6D=True)

    pysxt_dict = {
            'line_bb_dipole_not_cancelled': pysxt_line,
            'line_bb_dipole_cancelled': pysxt_line_bb_dipole_cancelled,
            'particle_on_closed_orbit': part_on_CO}

    if pickle_lines_in_folder is not None:
        pysix_fol_name = pickle_lines_in_folder
        os.makedirs(pysix_fol_name, exist_ok=True)

        with open(pysix_fol_name + "/line_bb_dipole_not_cancelled.pkl", "wb") as fid:
            pickle.dump(pysxt_line.to_dict(keepextra=True), fid)

        with open(pysix_fol_name + "/line_bb_dipole_cancelled.pkl", "wb") as fid:
            pickle.dump(pysxt_line_bb_dipole_cancelled.to_dict(keepextra=True), fid)

        with open(pysix_fol_name + "/particle_on_closed_orbit.pkl", "wb") as fid:
            pickle.dump(part_on_CO.to_dict(), fid)

    return pysxt_dict



import os

def make_links(links_dict, force=False):
    for kk in links_dict.keys():
        if force:
            if os.path.exists(kk):
                os.remove(kk)
        os.symlink(links_dict[kk], kk)

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
                mad_b4.input(f'const {nn}={b2_const[nn]}')

    # %% INDEPENDENT
    b2_indep=var_dicts_b2['independent_variables']
    b4_indep=var_dicts_b4['independent_variables']
    for nn in b2_indep.keys():
        mad_b4.input(f'{nn}={b2_indep[nn]}')

    # %% DEPENDENT
    b2_dep=var_dicts_b2['dependent_variables_expr']
    b4_dep=var_dicts_b4['dependent_variables_expr']
    for nn in b2_dep.keys():
        mad_b4.input(f'{nn}:={str(b2_dep[nn])}')

    # bv_aux and my my lhcbeam need to be defined explicitly
    mad_b4.input(f'bv_aux=-1')
    mad_b4.input(f'mylhcbeam=4')

    # Attach beam
    mad_b4.input(str(mad_b2.sequence['lhcb2'].beam))
    mad_b4.use('lhcb2')
    mad_b4.sequence['lhcb2'].beam['bv']=1

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



def checks_on_parameter_dict(params):

    assert params['par_nco_IP5']==params['par_nco_IP1']
    assert 'par_beam_norm_emit' in params
    print('Checks on paramter dict passed!')

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

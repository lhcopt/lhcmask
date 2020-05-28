import os

def make_links(links_dict, force=False):
    for kk in links_dict.keys():
        if force:
            if os.path.exists(kk):
                os.remove(kk)
        os.symlink(links_dict[kk], kk)

def checks_on_parameter_dict(params):

    assert params['par_nco_IP5']==params['par_nco_IP1']
    assert 'par_beam_norm_emit' in params
    print('Checks on paramter dict passed!')

def check_twiss_value(twiss_df, element_name, keyword, target, tol):
    assert abs(twiss_df.loc[element_name][keyword] - target) < tol,\
                f'Check not passes on {keyword} at {element_name}'

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
                f'Check not passes on {plane} separation at {element_name}'

def check_separations_against_madvars(checks, twiss_df_b1, twiss_df_b2, variables_dict):
    for cc in checks:
        tol = cc['tol']
        target = variables_dict['all_variables_val'][cc['varname']]
        check_separation_value(twiss_df_b1, twiss_df_b2, cc['element_name'],
                cc['plane'], target, cc['tol'])

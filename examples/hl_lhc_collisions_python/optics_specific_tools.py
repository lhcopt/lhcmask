import pymasktools as pmt

def check_beta_at_ips_against_madvars(beam, twiss_df, variable_dicts, tol):
    twiss_value_checks=[]
    for ip in [1,2,5,8]:
        for plane in ['x', 'y']:
            twiss_value_checks.append({
                    'element_name': f'ip{ip}:1',
                    'keyword': f'bet{plane}',
                    'varname': f'bet{plane}ip{ip}b{beam}',
                    'tol': 1e-3})

    pmt.check_twiss_against_madvars(twiss_value_checks, twiss_df, variable_dicts)

def check_separations_at_ips_against_madvars(twiss_df_b1, twiss_df_b2,
        variables_dict, tol):

    separations_to_check = []
    for ip in [1,2,5,8]:
        for plane in ['x', 'y']:
            separations_to_check.append({
                    'element_name': f'ip{ip}:1',
                    'plane': plane,
                    # knobs like on_sep1h, onsep8v etc
                    'varname': f'on_sep{ip}'+{'x':'h', 'y':'v'}[plane],
                    'tol': tol})
    pmt.check_separations_against_madvars(separations_to_check,
            twiss_df_b1, twiss_df_b2, variables_dict)

def twiss_and_check(mad, sequences_to_check, twiss_fname,
        check_betas_at_ips=True, check_separations_at_ips=True):
    var_dict = mad.get_variables_dicts()
    twiss_dfs = {}
    for ss in sequences_to_check:
        mad.use(ss)
        mad.twiss()
        tdf = mad.get_table_df('twiss')
        twiss_dfs[ss] = tdf

    for ss in sequences_to_check:
        tt = twiss_dfs[ss]
        if twiss_fname is not None:
            tt.to_parquet(twiss_fname + f'_{ss}.parquet')

    if check_betas_at_ips:
        for ss in sequences_to_check:
            tt = twiss_dfs[ss]
            check_beta_at_ips_against_madvars(beam=ss[-1],
                    twiss_df=tt,
                    variable_dicts=var_dict,
                    tol=1e-3)
        print('IP beta test against knobs passed!')

    if check_separations_at_ips:
        twiss_df_b1 = twiss_dfs['lhcb1']
        twiss_df_b2 = twiss_dfs['lhcb2']
        check_separations_at_ips_against_madvars(twiss_df_b1, twiss_df_b2,
                var_dict, tol=1e-6)
        print('IP separation test against knobs passed!')

    return twiss_dfs, var_dict

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


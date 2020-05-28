import pymasktools as pmt

def check_beta_at_ips_against_knobs(beam, twiss_df, variable_dicts, tol):
    twiss_value_checks=[]
    for ip in [1,2,5,8]:
        for plane in ['x', 'y']:
            twiss_value_checks.append({
                    'element_name': f'ip{ip}:1',
                    'keyword': f'bet{plane}',
                    'knob': f'bet{plane}ip{ip}b{beam}',
                    'tol': 1e-3})

    pmt.check_twiss_against_knobs(twiss_value_checks, twiss_df, variable_dicts)

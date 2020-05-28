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

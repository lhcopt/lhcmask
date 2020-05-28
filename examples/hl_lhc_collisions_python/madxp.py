from cpymad.madx import Madx

class Madxp(Madx):
    pass

    def set_variables_from_dict(self, params):
        for nn in params.keys():
            self.input(f'{nn}={params[nn]};')

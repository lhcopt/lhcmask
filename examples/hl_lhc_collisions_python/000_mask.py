import os

#from cpymad.madx import Madx
from madxp import Madxp as Madx

os.system('rm tracking_tools errors')

os.symlink('/afs/cern.ch/eng/tracking-tools', 'tracking_tools')
os.symlink('tracking_tools/errors', 'errors')

mad = Madx()

mad.input('mylhcbeam = 1')
mad.call('hl14_thin.madx')
mad.call('hl14_collision_optics.madx')

from parameters import parameters
mad.set_variables_from_dict(params=parameters)


import os

#from cpymad.madx import Madx
from madxp import Madxp as Madx
from tools import make_links

make_links(force=True, links_dict={
    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
    'modules': 'tracking_tools/modules',
    'tools': 'tracking_tools/tools',
    'beambeam_macros': 'tracking_tools/beambeam_macros',
    'errors': 'tracking_tools/errors'})


mad = Madx()

mad.input('mylhcbeam = 1')
mad.call('hl14_thin.madx')
mad.call('hl14_collision_optics.madx')

from parameters import parameters
mad.set_variables_from_dict(params=parameters)


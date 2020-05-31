import os

#from cpymad.madx import Madx
from madxp import Madxp as Madx
from pymasktools import make_links

make_links(force=True, links_dict={
    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
    'modules': 'tracking_tools/modules',
    'tools': 'tracking_tools/tools',
    'beambeam_macros': 'tracking_tools/beambeam_macros',
    'errors': 'tracking_tools/errors'})


mad = Madx()

# Build sequence
mad.input('mylhcbeam = 4')
mad.call('hl14_thin.madx')

# Make optics
mad.call('hl14_collision_optics.madx')

# Load parameters
from parameters import parameters
parameters['par_on_bb_switch'] = 0
mad.set_variables_from_dict(params=parameters)

# Call mask modules
mad.call("modules/module_01_orbit.madx")
mad.call("modules/module_02_lumilevel.madx")
mad.call("module_03_beambeam.madx")
# mad.use('lhcb2')
# mad.input('on_disp=0.')
# mad.use('lhcb2')
# mad.twiss();
# mad.use('lhcb2')


mad.call("modules/module_04_errors.madx")
mad.call("modules/module_05_tuning.madx")
mad.call("modules/module_06_generate.madx")

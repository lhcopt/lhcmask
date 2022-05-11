 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from cpymad.madx import Madx

from pymask import pymasktools as pmt 


# Importing LHC sequences:
mad = Madx()
mad.option(echo = True, warn = True)
mad.call(file='/afs/cern.ch/eng/lhc/optics/runII/2018/lhc_as-built.seq')
mad.call(file='/afs/cern.ch/eng/lhc/optics/runIII/RunIII_dev/2021_V6/PROTON/opticsfile.30')

 # Adding beam and extracting reference twiss table
mad.command.beam(particle='proton',charge=1,npart=1,energy=7000)
mad.sequence.lhcb1.use()
mad.twiss()

twiss_ref = mad.table['twiss'].dframe()
arc_ref   = twiss_ref.loc['s.arc.12.b1':'e.arc.12.b1']

#==================
# Removing sextupoles
#==================
all_sextupoles = arc_ref.groupby('keyword').get_group('sextupole')
testStr1 = pmt.seqedit(mad,seq_name='lhcb1',
                        editing = { 'mode'      : 'remove',
                                    'element'   : list(all_sextupoles.index)})

mad.twiss()
arc_no_sext = mad.table['twiss'].dframe().loc['s.arc.12.b1':'e.arc.12.b1']



#==================
# Installing markers
#==================
testStr2 = pmt.seqedit(mad,seq_name='lhcb1',
                        editing = { 'mode'      : 'install',
                                    'element'   : [f'marker.{sext}' for sext in all_sextupoles.index],
                                    'class'     : 'marker',
                                    'at'        : [mad.elements[sext]['at'] for sext in all_sextupoles.index],
                                    'from'      : [mad.elements[sext]['from'] for sext in all_sextupoles.index]})

mad.twiss()
arc_markers = mad.table['twiss'].dframe().loc['s.arc.12.b1':'e.arc.12.b1']


#==================
# Replacing with original sextupoles
#==================
testStr3 = pmt.seqedit(mad,seq_name='lhcb1',
                        editing = { 'mode'      : 'replace',
                                    'element'   : [f'marker.{sext}' for sext in all_sextupoles.index],
                                    'by'        : list(all_sextupoles.index)})

mad.twiss()
arc_recovered = mad.table['twiss'].dframe().loc['s.arc.12.b1':'e.arc.12.b1']



#==================
# Testing that everything went smoothly
#==================
# These should be different
assert(~arc_ref.equals(arc_no_sext))

# These should be identical
assert(arc_ref.equals(arc_recovered))
import sys
import numpy as np
import xtrack as xt
import pandas as pd
sys.path.append('./modules')
import pymask as pm

Madx = pm.Madxp
mad = Madx(command_log="mad_final.log")
mad.call("final_seq.madx")
mad.use(sequence="lhcb1")
mad.twiss()
mad.readtable(file="final_errors.tfs", table="errtab")
mad.seterr(table="errtab")
mad.set(format=".15g")
mad.twiss(rmatrix = True)

sequence_to_track = 'lhcb1'
bb_df_track = pd.read_pickle('bb_df_track.pickle')
optics_and_co_at_start_ring_from_madx = pm.get_optics_and_orbit_at_start_ring(
    mad, sequence_to_track, skip_mad_use=True)
#with open('./optics_orbit_at_start_ring_from_madx.json', 'w') as fid:
#    json.dump(optics_and_co_at_start_ring_from_madx, fid, cls=pm.JEncoder

########################                                                           
# Generate xtrack line #                                                           
########################                                                           
#if enable_bb_legacy:
#    print('xtrack line is not generated with bb legacy macros')
#    else:

tracker, line_bb_for_tracking_dict = pm.generate_xsuite_line(mad,
                        sequence_to_track, bb_df_track,
                        optics_and_co_at_start_ring_from_madx,
                        folder_name = './xsuite_lines',
                        skip_mad_use=True,
                        prepare_line_for_xtrack=True)

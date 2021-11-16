# Requires dataframes to be saved manually after executing 000_pymask.py

import sys
import json

import numpy as np

import xtrack as xt
from scipy.optimize import fsolve

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

separation_given_wrt_closed_orbit_4D = True
tw = mad.table['twiss']

seq_name = 'lhcb1'
mad_beam =  mad.sequence[seq_name].beam

import pandas as pd
bb_df_track = pd.read_pickle('./bb_df_track.pkl')
with open('./optics_orbit_at_start_ring_from_madx.json', 'r') as fid:
    optics_and_co_at_start_ring_from_madx = json.load(fid)

pm.generate_xsuite_line(mad, seq_name, bb_df_track,
                    optics_and_co_at_start_ring_from_madx,
                    folder_name = 'xsuite_lines',
                    skip_mad_use=False, prepare_line_for_xtrack=True)


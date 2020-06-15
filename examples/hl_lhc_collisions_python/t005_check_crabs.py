import shutil
import pickle
import numpy as np

import pysixtrack
import sixtracktools

# Test b1 
path_test = './'
type_test = 'sixtrack'

def prepare_line(path, input_type):

    if input_type == 'pysixtrack':
        # Load pysixtrack machine 
        with open(path, "rb") as fid:
            ltest = pysixtrack.Line.from_dict(pickle.load(fid))
    elif input_type == 'sixtrack':
        print('Build pysixtrack from sixtrack input:')
        sixinput_test = sixtracktools.sixinput.SixInput(path)
        ltest = pysixtrack.Line.from_sixinput(sixinput_test)
        print('Done')
    else:
        raise ValueError('What?!')

    return ltest

# Load
ltest = prepare_line(path_test, type_test)

# bb6d_elems, bb6d_names = ltest.get_elements_of_type(pysixtrack.elements.BeamBeam6D)
# 
ip_choice = 5
# 
# bb6d_elems_ip, bb6d_names_ip = list(zip(*[(ee, nn) for ee, nn in zip(bb6d_elems, bb6d_names)
#                                   if f'{ip_choice}b' in nn]))i

s_all_elems = ltest.get_s_elements()
names_all_elems = ltest.element_names
all_elems = ltest.elements

bb_names = []
bb_s = []
bb_elems = []
for nn, ss, ee in zip(names_all_elems, s_all_elems, all_elems):
    if nn.startswith('bb_ho'):
        if f'{ip_choice}b' in nn:
            bb_names.append(nn)
            bb_s.append(ss)
            bb_elems.append(ee)

iho = np.where(['.c' in nn for nn in bb_names])[0][0]
sho = bb_s[iho]
s_rel = np.array(bb_s) - sho

import matplotlib.pyplot as plt
plt.close('all')
plt.figure()
plt.plot(s_rel, [bb.y_bb_co for bb in bb_elems], 'o')

L_lhc = 27e3
h_cc = 35640
phi = 250e-6
phi_c = -190e-6

# For B1 as weak beam
X_crab = -phi_c * L_lhc / (2*np.pi*h_cc) *np.sin(2*np.pi*h_cc/L_lhc*2*s_rel)
X_no_crab = -phi * s_rel
X1_orbit = phi * s_rel
plt.plot(s_rel, X_no_crab + X_crab - X1_orbit, '*')
#plt.figure()
#plt.plot(s_rel, np.array([bb.y_bb_co for bb in bb_elems])
#                    /(2*s_rel), 'o')




plt.show()

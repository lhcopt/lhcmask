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
bb_index = []
for ii, (nn, ss, ee) in enumerate(zip(names_all_elems, s_all_elems, all_elems)):
    if nn.startswith('bb_ho'):
        if f'{ip_choice}b' in nn:
            bb_names.append(nn)
            bb_s.append(ss)
            bb_elems.append(ee)
            bb_index.append(ii)

iho = np.where(['.c' in nn for nn in bb_names])[0][0]
sho = bb_s[iho]
s_rel = np.array(bb_s) - sho

import matplotlib.pyplot as plt
plt.close('all')
fig1 = plt.figure(1)
axcrab = fig1.add_subplot(111)

L_lhc = 27e3
h_cc = 35640
phi = 250e-6
phi_c = -190e-6

# For B1 as weak beam
Y_crab = -phi_c * L_lhc / (2*np.pi*h_cc) *np.sin(2*np.pi*h_cc/L_lhc*2*s_rel)
Y_no_crab = -phi * s_rel
Y1_orbit = phi * s_rel
axcrab.plot(s_rel, np.array([bb.y_bb_co for bb in bb_elems])+Y1_orbit, 'o')
plt.plot(s_rel, Y_no_crab + Y_crab, '*')
#plt.figure()
#plt.plot(s_rel, np.array([bb.y_bb_co for bb in bb_elems])
#                    /(2*s_rel), 'o')


# Chek crabs in weak beam
import pickle
with open('./optics_orbit_at_start_ring.pkl', 'rb') as fid:
    ddd = pickle.load(fid)

ddd['p0c'] =  ddd['p0c_eV']
partco = pysixtrack.Particles.from_dict(ddd)


# Switch off all beam-beam lenses
bb_all, _ = ltest.get_elements_of_type([pysixtrack.elements.BeamBeam4D,
                                            pysixtrack.elements.BeamBeam6D])
for bb in bb_all: bb.enabled = False

# # Switch off all beam-beam lenses
crabs, _ = ltest.get_elements_of_type([pysixtrack.elements.RFMultipole])
#for cc in crabs:
#    cc.pn = [-90]
#    cc.ps = [-90]

# for cc in crabs:
#     cc.ksl = [0]
#     cc.knl = [0]

z_slices = s_rel * 2.0
partco.zeta += z_slices
partco.x += 0*z_slices
partco.s += 0*z_slices
partco.y += 0*z_slices
partco.px += 0*z_slices
partco.py += 0*z_slices
partco.delta += 0*z_slices


list_co = ltest.track_elem_by_elem(partco)
list_co2 = ltest.track_elem_by_elem(partco)

plt.figure(2)
axcox = plt.subplot(2,1,1)
axcoy = plt.subplot(2,1,2, sharex=axcox)
plt.suptitle('Check closed orbit 2 turns')
axcox.plot([pp.s[iho] for pp in list_co],  [pp.x[iho] for pp in list_co])
axcox.plot([pp.s[iho] for pp in list_co],  [pp.x[iho] for pp in list_co2])

axcoy.plot([pp.s[iho] for pp in list_co],  [pp.y[iho] for pp in list_co])
axcoy.plot([pp.s[iho] for pp in list_co],  [pp.y[iho] for pp in list_co2])

plt.figure(20)
axcopx = plt.subplot(2,1,1, sharex=axcox)
axcopy = plt.subplot(2,1,2, sharex=axcox)
plt.suptitle('Check closed orbit 2 turns')
axcopx.plot([pp.s[iho] for pp in list_co],  [pp.px[iho] for pp in list_co])
axcopx.plot([pp.s[iho] for pp in list_co],  [pp.px[iho] for pp in list_co2])

axcopy.plot([pp.s[iho] for pp in list_co],  [pp.py[iho] for pp in list_co])
axcopy.plot([pp.s[iho] for pp in list_co],  [pp.py[iho] for pp in list_co2])

plt.figure(3)
axcox = plt.subplot(2,1,1)
axcoy = plt.subplot(2,1,2, sharex=axcox)
plt.suptitle('Check crab')
for iz, zz in enumerate(z_slices):
    axcox.plot([pp.s[iz] for pp in list_co],  [pp.x[iz] for pp in list_co])
    axcoy.plot([pp.s[iz] for pp in list_co],  [pp.y[iz] for pp in list_co])

y_lenses = []
for ibb, bb in enumerate(bb_elems):
    y_lenses.append(list_co[bb_index[ibb]].y[ibb])

axcrab.plot(s_rel, y_lenses, 'o')

# part_list = []
# for nn, ee in zip(ltest.element_names, ltest.elements):
#     print(nn)
#     ee.track(partco)
#     part_list.append(partco.copy())
#     print(partco)
#     if np.any(np.abs(partco.delta)> 1e-3):
#         prrrr
#     if np.any(np.abs(partco.zeta)> 1.):
#         prrrr
#     if np.isnan(np.sum(partco.x)):
#         prrr
#     if np.isnan(np.sum(partco.px)):
#         prrr
#     if np.isnan(np.sum(partco.y)):
#         prrr
#     if np.isnan(np.sum(partco.py)):
#         prrr

plt.show()

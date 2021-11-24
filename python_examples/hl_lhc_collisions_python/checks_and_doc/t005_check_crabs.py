import shutil
import numpy as np
import json

import xtrack as xt
import xpart as xp
import xfields as xf
import sixtracktools

L_lhc = 27e3
h_cc = 35640

# Mad-X knobs
Phi = 250e-6
Phi_c = -190e-6

# B1 ip5
beam_track = 'b1'
ip_choice = 5
plane = 'y'
phi_weak = Phi
phi_c_weak = Phi_c

# # B1 ip1
# beam_track = 'b1'
# ip_choice = 1
# plane = 'x'
# phi_weak = Phi
# phi_c_weak = Phi_c

# # B4 ip5
# beam_track = 'b4'
# ip_choice = 5
# plane = 'y'
# phi_weak = Phi
# phi_c_weak = Phi_c

# # B4 ip1
# beam_track = 'b4'
# ip_choice = 1
# plane = 'x'
# phi_weak = -Phi
# phi_c_weak = -Phi_c

phi_strong = -phi_weak
phi_c_strong = -phi_c_weak

path_test = '../'
type_test = 'sixtrack'

#path_test = '../xsuite_lines/line_bb_dipole_not_cancelled.json'
#type_test = 'xsuite'

def prepare_line(path, input_type):

    if input_type == 'xsuite':
        # Load xsuite machine 
        with open(path, 'r') as fid:
            ltest = xt.Line.from_dict(json.load(fid))
    elif input_type == 'sixtrack':
        print('Build xsuite from sixtrack input:')
        sixinput_test = sixtracktools.sixinput.SixInput(path)
        ltest = xt.Line.from_sixinput(sixinput_test)
        print('Done')
    else:
        raise ValueError('What?!')

    return ltest

# Load
ltest = prepare_line(path_test, type_test)


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


# For B1 as weak beam (and we plot B2 as stron beam)
# phi_c_strong = dR_crab / dz
# phi_strong = dR_crab / ds
R_crab = phi_c_strong * L_lhc / (2*np.pi*h_cc) *np.sin(2*np.pi*h_cc/L_lhc*2*s_rel)
R_no_crab = phi_strong * s_rel
R1_orbit = phi_weak * s_rel

axcrab.plot(s_rel, np.array(
    [getattr(bb, f'delta_{plane}') for bb in bb_elems])+R1_orbit,
    'o', color = 'r', alpha=.5, label='strong xsuite')
plt.plot(s_rel, R_no_crab + R_crab, '*', color='darkred', label='strong formula')

#axcrab.plot(s_rel, np.array([bb.y_bb_co for bb in bb_elems]), 'o', color='r', alpha=.5)
#plt.plot(s_rel, Y_no_crab + Y_crab - Y1_orbit, '*', color='darkred')

#plt.plot(s_rel, Y_no_crab, 'xr')
#plt.plot(s_rel, Y_crab, 'xb')
#plt.figure()
#plt.plot(s_rel, np.array([bb.y_bb_co for bb in bb_elems])
#                    /(2*s_rel), 'o')


# Chek crabs in weak beam

# Switch off all beam-beam lenses
bb_all, _ = ltest.get_elements_of_type([xf.BeamBeamBiGaussian2D,
                                        xf.BeamBeamBiGaussian3D])
for bb in bb_all:
    bb.q0 = 0

# # Switch off all beam-beam lenses
crabs, crab_names = ltest.get_elements_of_type([xt.RFMultipole])
#for cc in crabs:
#    cc.pn = [-90]
#    cc.ps = [-90]

# for cc in crabs:
#     cc.knl[0] *= -1

with open('../optics_orbit_at_start_ring_from_madx.json', 'r') as fid:
    ddd = json.load(fid)
partco = xp.Particles.from_dict(ddd['particle_on_madx_co'])
z_slices = s_rel * 2.0
partco = xp.build_particles(particle_ref=partco, zeta=z_slices)


list_co = ltest.slow_track_elem_by_elem(partco)
list_co2 = ltest.slow_track_elem_by_elem(partco)

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

# For each s_lens, we find the transverse position of the weak beam slice 
# that collides with the sycnhronous particle of the strong
r_lenses = []
for ibb, bb in enumerate(bb_elems):
    r_lenses.append(getattr(list_co[bb_index[ibb]], plane)[ibb])

axcrab.plot(s_rel, r_lenses, 'o', color='b', alpha=.5, label= 'weak xsuite')
Rw_crab = phi_c_weak * L_lhc / (2*np.pi*h_cc) *np.sin(2*np.pi*h_cc/L_lhc*2*s_rel)
Rw_no_crab = phi_weak * s_rel
axcrab.plot(s_rel, Rw_crab + Rw_no_crab, '*', color='darkblue',
        label='weak formula')
axcrab.legend(loc='best')

# Check crab bump
import pandas as pd
z_crab_twiss = 0.075

if beam_track == 'b1':
    crab_df = pd.read_parquet(f'../twiss_z_crab_{z_crab_twiss:.5f}_seq_lhcb1.parquet')
    s_twiss = crab_df.s.values
    x_twiss = crab_df.x.values
    y_twiss = crab_df.y.values
    px_twiss = crab_df.px.values
    py_twiss = crab_df.py.values
    z_crab_track = z_crab_twiss
elif beam_track == 'b4':
    crab_df = pd.read_parquet(f'../twiss_z_crab_{z_crab_twiss:.5f}_seq_lhcb2.parquet')
    s_twiss = -crab_df.s.values[::-1]
    s_twiss -= s_twiss[0]
    x_twiss = -crab_df.x.values[::-1]
    y_twiss = crab_df.y.values[::-1]
    px_twiss = crab_df.px.values[::-1]
    py_twiss = -crab_df.py.values[::-1]
    z_crab_track = -z_crab_twiss
else:
    raise ValueError('????!!!!!')

figcb = plt.figure(4)
axcbx = figcb.add_subplot(2,1,1)
axcby = figcb.add_subplot(2,1,2, sharex=axcbx)

axcbx.plot(s_twiss, x_twiss)
axcby.plot(s_twiss, y_twiss)

part = xp.Particles.from_dict(ddd['particle_on_madx_co'])

z_test = np.array([0, z_crab_track])
part = xp.build_particles(particle_ref=part,
        zeta = z_test,
        x = 0*z_test + np.array([0, x_twiss[0]]),
        y = 0*z_test + np.array([0, y_twiss[0]]),
        px = 0*z_test + np.array([0, px_twiss[0]]),
        py = 0*z_test + np.array([0, py_twiss[0]]),
        delta = 0*z_test)
list_track = ltest.slow_track_elem_by_elem(part)

axcbx.plot([pp.s[0] for pp in list_track],  [pp.x[1] - pp.x[0] for pp in list_track])
axcby.plot([pp.s[0] for pp in list_track],  [pp.y[1] - pp.y[0] for pp in list_track])

plt.show()

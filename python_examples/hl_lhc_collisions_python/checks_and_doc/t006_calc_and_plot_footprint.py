import pickle
import numpy as np
import NAFFlib
import helpers as hp
import footprint
import matplotlib.pyplot as plt

import pysixtrack
import sixtracktools

track_with = 'PySixtrack'
track_with = 'Sixtracklib'
#device = 'opencl:1.0'
device = None

opt_start_ring_fname = '../optics_orbit_at_start_ring.pkl'

epsn_x = 2.5e-6
epsn_y = 2.5e-6
r_max_sigma = 5.
N_r_footp = 20
N_theta_footp = 10

n_turns = 100

def prepare_line(path, input_type):

    if input_type == 'pysixtrack':
        # Load pysixtrack machine 
        with open(path, "rb") as fid:
            ltest = pysixtrack.Line.from_dict(
                    pickle.load(fid))
    elif input_type == 'sixtrack':
        print('Build pysixtrack from sixtrack input:')
        sixinput_test = sixtracktools.sixinput.SixInput(path)
        ltest = pysixtrack.Line.from_sixinput(sixinput_test)
        print('Done')
    else:
        raise ValueError('What?!')

    return ltest

line = prepare_line('../', input_type='sixtrack')


with open('../optics_orbit_at_start_ring.pkl', 'rb') as fid:
    ddd = pickle.load(fid)
ddd['p0c'] =  ddd['p0c_eV']

partCO = pysixtrack.Particles.from_dict(ddd)

# line.disable_beambeam()
part = partCO.copy()

line.beambeam_store_closed_orbit_and_dipolar_kicks(
        partCO.copy(),
        separation_given_wrt_closed_orbit_4D=True,
        separation_given_wrt_closed_orbit_6D=True)



beta_x = ddd['betx']
beta_y = ddd['bety']

sigmax = np.sqrt(beta_x * epsn_x / part.beta0 / part.gamma0)
sigmay = np.sqrt(beta_y * epsn_y / part.beta0 / part.gamma0)

xy_norm = footprint.initial_xy_polar(r_min=1e-2, r_max=r_max_sigma, r_N=N_r_footp + 1,
                                     theta_min=np.pi / 100, theta_max=np.pi / 2 - np.pi / 100,
                                     theta_N=N_theta_footp)

DpxDpy_wrt_CO = np.zeros_like(xy_norm)

for ii in range(xy_norm.shape[0]):
    for jj in range(xy_norm.shape[1]):

        DpxDpy_wrt_CO[ii, jj, 0] = xy_norm[ii, jj, 0] * np.sqrt(epsn_x / part.beta0 / part.gamma0 / beta_x)
        DpxDpy_wrt_CO[ii, jj, 1] = xy_norm[ii, jj, 1] * np.sqrt(epsn_y / part.beta0 / part.gamma0 / beta_y)


if track_with == 'PySixtrack':

    part = partCO.copy()

    x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_pysixtrack(
        line, part=part, Dx_wrt_CO_m=0., Dpx_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 0].flatten(),
        Dy_wrt_CO_m=0, Dpy_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 1].flatten(),
        Dsigma_wrt_CO_m=0., Ddelta_wrt_CO=0., n_turns=n_turns, verbose=True)

    info = track_with

elif track_with == 'Sixtrack':
    x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(
        partCO=partCO, Dx_wrt_CO_m=0., Dpx_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 0].flatten(),
        Dy_wrt_CO_m=0, Dpy_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 1].flatten(),
        Dsigma_wrt_CO_m=0., Ddelta_wrt_CO=0., n_turns=n_turns,
        input_folder=sixtrack_input_folder)
    info = track_with

elif track_with == 'Sixtracklib':
    x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtracklib(
        line=line, partCO=partCO, Dx_wrt_CO_m=0., Dpx_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 0].flatten(),
        Dy_wrt_CO_m=0., Dpy_wrt_CO_rad=DpxDpy_wrt_CO[:, :, 1].flatten(),
        Dsigma_wrt_CO_m=0., Ddelta_wrt_CO=0., n_turns=n_turns, device=device)
    info = track_with
    if device is None:
        info += ' (CPU)'
    else:
        info += ' (GPU %s)'%device
else:
    raise ValueError('What?!')

n_part = x_tbt.shape[1]
Qx = np.zeros(n_part)
Qy = np.zeros(n_part)

for i_part in range(n_part):
    Qx[i_part] = NAFFlib.get_tune(x_tbt[:, i_part])
    Qy[i_part] = NAFFlib.get_tune(y_tbt[:, i_part])

Qxy_fp = np.zeros_like(xy_norm)

Qxy_fp[:, :, 0] = np.reshape(Qx, Qxy_fp[:, :, 0].shape)
Qxy_fp[:, :, 1] = np.reshape(Qy, Qxy_fp[:, :, 1].shape)

plt.close('all')

fig3 = plt.figure(3)
axcoord = fig3.add_subplot(1, 1, 1)
footprint.draw_footprint(xy_norm, axis_object=axcoord, linewidth = 1)
axcoord.set_xlim(right=np.max(xy_norm[:, :, 0]))
axcoord.set_ylim(top=np.max(xy_norm[:, :, 1]))

fig4 = plt.figure(4)
axFP = fig4.add_subplot(1, 1, 1)
footprint.draw_footprint(Qxy_fp, axis_object=axFP, linewidth = 1)
# axFP.set_xlim(right=np.max(Qxy_fp[:, :, 0]))
# axFP.set_ylim(top=np.max(Qxy_fp[:, :, 1]))
fig4.suptitle(info)
plt.show()

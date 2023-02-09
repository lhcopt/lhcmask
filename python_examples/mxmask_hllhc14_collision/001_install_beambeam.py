import json
import numpy as np

import xtrack as xt
import xobjects as xo

with open('collider_00_from_mad.json', 'r') as fid:
    collider = xt.Multiline.from_dict(json.load(fid))

collider.install_beambeam_interactions(
    clockwise_line='lhcb1',
    anticlockwise_line='lhcb2',
    ip_names=['ip1', 'ip2', 'ip5', 'ip8'],
    num_long_range_encounters_per_side=[25, 20, 25, 20],
    num_slices_head_on=11,
    harmonic_number=35640,
    bunch_spacing_buckets=10,
    sigmaz=0.076)

collider.build_trackers()

with open('collider_01_bb_off.json', 'w') as fid:
    dct = collider.to_dict()
    json.dump(dct, fid, cls=xo.JEncoder)


###### Checks ######

collider_before_save = collider
dct = collider.to_dict()
collider = xt.Multiline.from_dict(dct)
collider.build_trackers()

assert collider._bb_config['dataframes']['clockwise'].shape == (
    collider_before_save._bb_config['dataframes']['clockwise'].shape)
assert collider._bb_config['dataframes']['anticlockwise'].shape == (
    collider_before_save._bb_config['dataframes']['anticlockwise'].shape)

assert (collider._bb_config['dataframes']['clockwise']['elementName'].iloc[50]
    == collider_before_save._bb_config['dataframes']['clockwise']['elementName'].iloc[50])
assert (collider._bb_config['dataframes']['anticlockwise']['elementName'].iloc[50]
    == collider_before_save._bb_config['dataframes']['anticlockwise']['elementName'].iloc[50])

# Put in some orbit
knobs = dict(on_x1=250, on_x5=-200, on_disp=1)

for kk, vv in knobs.items():
    collider.vars[kk] = vv

tw1_b1 = collider['lhcb1'].twiss(method='4d')
tw1_b2 = collider['lhcb2'].twiss(method='4d')

with open('collider_00_from_mad.json', 'r') as fid:
    collider_ref = xt.Multiline.from_dict(json.load(fid))

collider_ref.build_trackers()

for kk, vv in knobs.items():
    collider_ref.vars[kk] = vv

tw0_b1 = collider_ref['lhcb1'].twiss(method='4d')
tw0_b2 = collider_ref['lhcb2'].twiss(method='4d')

assert np.isclose(tw1_b1.qx, tw0_b1.qx, atol=1e-7, rtol=0)
assert np.isclose(tw1_b1.qy, tw0_b1.qy, atol=1e-7, rtol=0)
assert np.isclose(tw1_b2.qx, tw0_b2.qx, atol=1e-7, rtol=0)
assert np.isclose(tw1_b2.qy, tw0_b2.qy, atol=1e-7, rtol=0)

assert np.isclose(tw1_b1.dqx, tw0_b1.dqx, atol=1e-4, rtol=0)
assert np.isclose(tw1_b1.dqy, tw0_b1.dqy, atol=1e-4, rtol=0)
assert np.isclose(tw1_b2.dqx, tw0_b2.dqx, atol=1e-4, rtol=0)
assert np.isclose(tw1_b2.dqy, tw0_b2.dqy, atol=1e-4, rtol=0)

for ipn in [1, 2, 3, 4, 5, 6, 7, 8]:
    assert np.isclose(tw1_b1[f'ip{ipn}', 'betx'], tw0_b1[f'ip{ipn}', 'betx'], rtol=1e-5, atol=0)
    assert np.isclose(tw1_b1[f'ip{ipn}', 'bety'], tw0_b1[f'ip{ipn}', 'bety'], rtol=1e-5, atol=0)
    assert np.isclose(tw1_b2[f'ip{ipn}', 'betx'], tw0_b2[f'ip{ipn}', 'betx'], rtol=1e-5, atol=0)
    assert np.isclose(tw1_b2[f'ip{ipn}', 'bety'], tw0_b2[f'ip{ipn}', 'bety'], rtol=1e-5, atol=0)

    assert np.isclose(tw1_b1[f'ip{ipn}', 'px'], tw0_b1[f'ip{ipn}', 'px'], rtol=1e-9, atol=0)
    assert np.isclose(tw1_b1[f'ip{ipn}', 'py'], tw0_b1[f'ip{ipn}', 'py'], rtol=1e-9, atol=0)
    assert np.isclose(tw1_b2[f'ip{ipn}', 'px'], tw0_b2[f'ip{ipn}', 'px'], rtol=1e-9, atol=0)
    assert np.isclose(tw1_b2[f'ip{ipn}', 'py'], tw0_b2[f'ip{ipn}', 'py'], rtol=1e-9, atol=0)

    assert np.isclose(tw1_b1[f'ip{ipn}', 's'], tw0_b1[f'ip{ipn}', 's'], rtol=1e-10, atol=0)
    assert np.isclose(tw1_b2[f'ip{ipn}', 's'], tw0_b2[f'ip{ipn}', 's'], rtol=1e-10, atol=0)




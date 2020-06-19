import shutil
import pickle
import numpy as np

import pysixtrack
import sixtracktools

# Test b1 
path_test = './'
type_test = 'sixtrack'
path_ref = '../../examples/hl_lhc_collision'
type_ref = 'sixtrack'

# # Test b4 nobb sixtrack
# path_test = './'
# type_test = 'sixtrack'
# path_ref = '../hl_lhc_collision_nobb_b4'
# type_ref = 'sixtrack'
#
# # Test b4 nobb pysixtrack (not working for now)
# path_test = './pysixtrack/line_bb_dipole_not_cancelled.pkl'
# type_test = 'pysixtrack'
# path_ref = '../hl_lhc_collision_nobb_b4'
# type_ref = 'sixtrack'

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
lref = prepare_line(path_ref, type_ref)

original_length = ltest.get_length()
assert (lref.get_length() - original_length) < 1e-6

# Simplify the two machines
for ll in (ltest, lref):
    ll.remove_inactive_multipoles(inplace=True)
    ll.remove_zero_length_drifts(inplace=True)
    ll.merge_consecutive_drifts(inplace=True)

# Check that the two machines are identical
assert len(ltest) == len(lref)

assert (ltest.get_length() - original_length) < 1e-6
assert (lref.get_length() - original_length) < 1e-6


def norm(x):
    return np.sqrt(np.sum(np.array(x) ** 2))


for ii, (ee_test, ee_six, nn_test, nn_six) in enumerate(
    zip(ltest.elements, lref.elements, ltest.element_names, lref.element_names)
):
    assert type(ee_test) == type(ee_six)

    dtest = ee_test.to_dict(keepextra=True)
    dsix = ee_six.to_dict(keepextra=True)

    for kk in dtest.keys():

        # Check if they are identical
        if dtest[kk] == dsix[kk]:
            continue

        # Check if the relative error is small
        try:
            diff_rel = norm(np.array(dtest[kk]) - np.array(dsix[kk])) / norm(
                dtest[kk]
            )
        except ZeroDivisionError:
            diff_rel = 100.0
        if diff_rel < 3e-5:
            continue

        # Check if absolute error is small
        diff_abs = norm(np.array(dtest[kk]) - np.array(dsix[kk]))
        if diff_abs > 0:
            print(f"{nn_test}[{kk}] - test:{dtest[kk]} six:{dsix[kk]}")
        if diff_abs < 1e-12:
            continue

        # Exception: drift length (100 um tolerance)
        if isinstance(
            ee_test, (pysixtrack.elements.Drift, pysixtrack.elements.DriftExact)
        ):
            if kk == "length":
                if diff_abs < 1e-4:
                    continue

        # Exception: multipole lrad is not passed to sixtraxk
        if isinstance(ee_test, pysixtrack.elements.Multipole):
            if kk == "length":
                if np.abs(ee_test.hxl) + np.abs(ee_test.hyl) == 0.0:
                    continue

        # Exceptions BB4D (separations are recalculated)
        if isinstance(ee_test, pysixtrack.elements.BeamBeam4D):
            if kk == "x_bb":
                if diff_abs / dtest["sigma_x"] < 0.0001:
                    continue
            if kk == "y_bb":
                if diff_abs / dtest["sigma_y"] < 0.0001:
                    continue
            if kk == "sigma_x":
                if diff_rel < 1e-5:
                    continue
            if kk == "sigma_y":
                if diff_rel < 1e-5:
                    continue

        # Exceptions BB4D (angles and separations are recalculated)
        if isinstance(ee_test, pysixtrack.elements.BeamBeam6D):
            if kk == "alpha":
                if diff_abs < 10e-6:
                    continue
            if kk == "x_bb_co":
                if diff_abs / np.sqrt(dtest["sigma_11"]) < 0.004:
                    continue
            if kk == "y_bb_co":
                if diff_abs / np.sqrt(dtest["sigma_33"]) < 0.004:
                    continue

        # If it got here it means that no condition above is met
        raise ValueError("Too large discrepancy!")
print(
    """
*******************************************************************
 The line from test seq. and the line from reference are identical!
*******************************************************************
"""
)

import shutil
import pickle
import numpy as np

import xline
import sixtracktools

# Test b1 
path_test = '../'
type_test = 'sixtrack'
path_ref = '../../../examples/hl_lhc_collision'
type_ref = 'sixtrack'
rtol = 3e-5
atol = 1e-12
strict = False

# # Test b1 
# path_test = '../xline/line_bb_dipole_not_cancelled.json'
# type_test = 'xline'
# path_ref = '../'
# type_ref = 'sixtrack'
# rtol = 3e-7
# atol = 1e-100
# strict=True

# # Test b4 nobb sixtrack
# path_test = '../../'
# type_test = 'sixtrack'
# path_ref = '../../../examples/hl_lhc_collision_nobb_b4'
# type_ref = 'sixtrack'

# # Test b4 nobb xline (not working for now)
# path_test = './xline/line_bb_dipole_not_cancelled.pkl'
# type_test = 'xline'
# path_ref = '../hl_lhc_collision_nobb_b4'
# type_ref = 'sixtrack'

def prepare_line(path, input_type):

    if input_type == 'xline':
        # Load xline machine 
        ltest = xline.Line.from_json(path_test)
    elif input_type == 'sixtrack':
        print('Build xline from sixtrack input:')
        sixinput_test = sixtracktools.sixinput.SixInput(path)
        ltest = xline.Line.from_sixinput(sixinput_test)
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
    ll.merge_consecutive_multipoles(inplace=True)

    # Remove inactive RFMultipoles and normalize phase
    for ii, ee in enumerate(ll.elements):
        if ee.__class__.__name__ == 'RFMultipole':
            if ( np.max(np.abs(ee.knl)) < 1e-20 and
                 np.max(np.abs(ee.ksl)) < 1e-20 and
                 abs(ee.voltage) < 1e-20):
                ll.element_names[ii] = None
                ll.elements[ii] = None
            # # Normalize phase
            # for kkll, pp in [[ee.knl, ee.pn],
            #                  [ee.ksl, ee.ps]]:
            #     for ii, vv in enumerate(kkll):
            #         if vv < 0:
            #             kkll[ii] = -vv
            #             pp[ii] += 180
            #         if pp[ii]>180:
            #             pp[ii] -= 360

    while None in ll.element_names:
        ll.element_names.remove(None)
    while None in ll.elements:
        ll.elements.remove(None)


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
    dref = ee_six.to_dict(keepextra=True)

    for kk in dtest.keys():

        # Check if they are identical
        if np.isscalar(dref[kk]) and dtest[kk] == dref[kk]:
            continue

        # Check if the relative error is small
        val_test = dtest[kk]
        val_ref = dref[kk]
        try:
            if not np.isscalar(val_ref) and len(val_ref) != len(val_test):
                    diff_rel = 100
                #lmin = min(len(val_ref), len(val_test))
                #val_test = val_test[:lmin]
                #val_ref = val_ref[:lmin]
            else:
                diff_rel = norm(np.array(val_test) - np.array(val_ref)) / norm(val_test)
        except ZeroDivisionError:
            diff_rel = 100.0
        if diff_rel < rtol:
            continue

        # Check if absolute error is small

        if not np.isscalar(val_ref) and len(val_ref) != len(val_test):
            diff_abs = 1000
        else:
            diff_abs = norm(np.array(val_test) - np.array(val_ref))
        if diff_abs > 0:
            print(f"{nn_test}[{kk}] - test:{dtest[kk]} six:{dref[kk]}")
        if diff_abs < atol:
            continue

        # Exception: drift length (100 um tolerance)
        if not(strict) and isinstance(
            ee_test, (xline.elements.Drift, xline.elements.DriftExact)
        ):
            if kk == "length":
                if diff_abs < 1e-4:
                    continue

        # Exception: multipole lrad is not passed to sixtraxk
        if isinstance(ee_test, xline.elements.Multipole):
            if kk == "length":
                if np.abs(ee_test.hxl) + np.abs(ee_test.hyl) == 0.0:
                    continue
            if kk == 'knl' or kk == 'ksl':
                val_ref = np.trim_zeros(val_ref)
                val_test= np.trim_zeros(val_test)
                if len(val_ref) != len(val_test):
                    lmin = min(len(val_ref), len(val_test))
                    if lmin < 10:
                        raise ValueError('Missing significant multipole strength')
                    else:
                        val_ref = val_ref[:lmin]
                        val_test = val_test[:lmin]
                if len(val_ref) == 0 and len(val_test) == 0:
                    continue
                else:
                    diff_rel = (norm(np.array(val_test) - np.array(val_ref))
                                /norm(val_test))
                    if diff_rel < rtol:
                        continue

        # Exception: correctors involved in lumi leveling
        passed_corr = False
        for nn_corr in [
            'mcbcv.5l8.b1', 'mcbyv.a4l8.b1', 'mcbxv.3l8',
            'mcbyv.4r8.b1', 'mcbyv.b5r8.b1',
            'mcbyh.b5l2.b1', 'mcbyh.4l2.b1', 'mcbxh.3l2', 'mcbyh.a4r2.b1',
            'mcbch.5r2.b1',
            'mcbcv.5l8.b2', 'mcbyv.a4l8.b2', 'mcbxv.3l8',
            'mcbyv.4r8.b2', 'mcbyv.b5r8.b2',
            'mcbyh.b5l2.b2', 'mcbyh.4l2.b2', 'mcbxh.3l2', 'mcbyh.a4r2.b2',
            'mcbch.5r2.b2', 'mcbch.a5r2.b2', 'mcbyh.4r2.b2', 'mcbxh.3r2',
            'mcbyh.a4l2.b2', 'mcbyh.5l2.b2', 'mcbyv.5r8.b2', 'mcbyv.a4r8.b2',
            'mcbxv.3r8', 'mcbyv.4l8.b2', 'mcbcv.b5l8.b2']:
            if nn_corr in nn_test and diff_rel < 1e-2:
                passed_corr = True
                break
        if not(strict) and  passed_corr:
            continue


        # Exceptions BB4D (separations are recalculated)
        if not(strict) and isinstance(ee_test, xline.elements.BeamBeam4D):
            if kk == "x_bb":
                if diff_abs / dtest["sigma_x"] < 0.01: # This is neede to accommodate different leveling routines (1% difference)
                    continue
            if kk == "y_bb":
                if diff_abs / dtest["sigma_y"] < 0.01:
                    continue
            if kk == "sigma_x":
                if diff_rel < 1e-5:
                    continue
            if kk == "sigma_y":
                if diff_rel < 1e-5:
                    continue

        # Exceptions BB4D (angles and separations are recalculated)
        if not(strict) and isinstance(ee_test, xline.elements.BeamBeam6D):
            if kk == "alpha":
                if diff_abs < 10e-6:
                    continue
            if kk == "x_bb_co":
                if diff_abs / np.sqrt(dtest["sigma_11"]) < 0.015:
                    continue
            if kk == "y_bb_co":
                if diff_abs / np.sqrt(dtest["sigma_33"]) < 0.015:
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

import time
import shutil
import pickle
import json
import numpy as np

import sixtracktools
import xtrack as xt
import xfields as xf




# Tests b1 with bb
tests = [
    {
        'test_name': 'B1 - pymask sixtrack input vs madx mask',
        'path_test': '../',
        'type_test': 'sixtrack',
        'path_ref': '../../../examples/hl_lhc_collision',
        'type_ref': 'sixtrack',
        'rtol': 3e-5,
        'atol': 1e-12,
        'strict': False,
    },
    {
        'test_name': 'B1 - pymask xsuite vs pymask sixtrack input',
        'path_test': '../xsuite_lines/line_bb_dipole_not_cancelled.json',
        'type_test': 'xsuite',
        'path_ref': '../',
        'type_ref': 'sixtrack',
        'rtol': 4e-7,
        'atol': 1e-100,
        'strict': True,
    }
]

# # Tests b4 no bb
# tests = [
#     {
#         'test_name': 'B4 - pymask sixtrack input vs madx mask',
#         'path_test': '../',
#         'type_test': 'sixtrack',
#         'path_ref': '../../../examples/hl_lhc_collision_nobb_b4',
#         'type_ref': 'sixtrack',
#         'rtol': 8e-5,
#         'atol': 1e-12,
#         'strict': False,
#     },
#     {
#         'test_name': 'B4 - pymask xsuite vs pymask sixtrack input',
#         'path_test': '../xsuite/line_bb_dipole_not_cancelled.json',
#         'type_test': 'xsuite',
#         'path_ref': '../',
#         'type_ref': 'sixtrack',
#         'rtol': 4e-7,
#         'atol': 1e-100,
#         'strict': True,
#     }
# ]

def norm(x):
    return np.sqrt(np.sum(np.array(x) ** 2))

def prepare_line(path, input_type):

    if input_type == 'xsuite':
        # Load machine 
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

for tt in tests:
    test_name = tt['test_name']
    path_test = tt['path_test']
    type_test = tt['type_test']
    path_ref = tt['path_ref']
    type_ref = tt['type_ref']
    rtol = tt['rtol']
    atol = tt['atol']
    strict = tt['strict']

    # Load
    ltest = prepare_line(path_test, type_test)
    lref = prepare_line(path_ref, type_ref)
    original_length = ltest.get_length()
    assert (lref.get_length() - original_length) < 1e-6

    # Simplify the two machines
    for ll in (ltest, lref):
        ll._var_management = None
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
                    ll.element_dict[ll.element_names[ii]] = None
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

    # Check that the two machines are identical
    assert len(ltest) == len(lref)

    assert (ltest.get_length() - original_length) < 1e-6
    assert (lref.get_length() - original_length) < 1e-6

    for ii, (ee_test, ee_six, nn_test, nn_six) in enumerate(
        zip(ltest.elements, lref.elements, ltest.element_names, lref.element_names)
    ):
        assert type(ee_test) == type(ee_six)

        dtest = ee_test.to_dict()
        dref = ee_six.to_dict()

        for kk in dtest.keys():

            # Check if they are identical
            if np.isscalar(dref[kk]) and dtest[kk] == dref[kk]:
                continue

            if isinstance(dref[kk], dict):
                if kk=='fieldmap':
                    continue
                if kk=='boost_parameters':
                    continue
                if kk=='Sigmas_0_star':
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
            if not(strict) and isinstance(ee_test, xt.Drift):
                if kk == "length":
                    if diff_abs < 1e-4:
                        continue

            # Exception: multipole lrad is not passed to sixtraxk
            if isinstance(ee_test, xt.Multipole):
                if kk == "length":
                    if np.abs(ee_test.hxl) + np.abs(ee_test.hyl) == 0.0:
                        continue
                if kk == "order":
                    # Checked through bal
                    continue
                if kk == 'knl' or kk == 'ksl' or kk == 'bal':
                    if len(val_ref) != len(val_test):
                        lmin = min(len(val_ref), len(val_test))
                        for vv in [val_ref,val_test]:
                            if len(vv)> lmin:
                                for oo in range(lmin, len(vv)):
                                    # we do not care about errors above 10
                                    if vv[oo] != 0 and oo < {'knl':10,
                                                         'ksl':10, 'bal':20}[kk]:
                                        raise ValueError(
                                            'Missing significant multipole strength')

                        val_ref = val_ref[:lmin]
                        val_test = val_test[:lmin]

                    if len(val_ref) == 0 and len(val_test) == 0:
                        continue
                    else:
                        diff_abs = norm(np.array(val_test) - np.array(val_ref))
                        diff_rel = diff_abs/norm(val_test)
                        if diff_rel < rtol:
                            continue
                        if diff_abs < atol:
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
            if not(strict) and isinstance(ee_test, xf.BeamBeamBiGaussian2D):
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
            if isinstance(ee_test, xf.BeamBeamBiGaussian2D):
                if kk == 'q0' or kk == 'n_particles':
                    # ambiguity due to old interface
                    if np.abs(ee_test.n_particles*ee_test.q0 -
                            ee_six.n_particles*ee_six.q0 ) < 1.: # charges
                        continue

            # Exceptions BB6D (angles and separations are recalculated)
            if not(strict) and isinstance(ee_test, xf.BeamBeamBiGaussian3D):
                if kk == "alpha":
                    if diff_abs < 10e-6:
                        continue
                if kk == "x_co" or kk == "x_bb_co" or kk == 'delta_x':
                    if diff_abs / np.sqrt(dtest["sigma_11"]) < 0.015:
                        continue
                if kk == "y_co" or kk == "y_bb_co" or kk == 'delta_y':
                    if diff_abs / np.sqrt(dtest["sigma_33"]) < 0.015:
                        continue
                if kk == "zeta_co":
                    if diff_abs <1e-5:
                        continue
                if kk == "delta_co":
                    if diff_abs <1e-5:
                        continue
                if kk == "px_co" or kk == 'py_co':
                    if diff_abs <30e-9:
                        continue

            # If it got here it means that no condition above is met
            raise ValueError("Too large discrepancy!")
    print(
       f"""
    *******************************************************************
     test: {test_name}

     The line from test seq. and the line from reference are identical!
    *******************************************************************
    """
    )

    time.sleep(1)

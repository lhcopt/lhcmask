import numpy as np

# Based on:
    # Fine tuning of coupling correction
    # S.Fartoukh March 2009
    # 2011/11/21 From SLHCV3.0 R. De Maria
    # 2016/04/27 R. De Maria & F.F Van der Veken
    #            work around mad-x bug on integer part of the tune in presence of coupling for ctap calculation

def coupling_measurement(mad, qx_integer, qy_integer,
        qx_fractional, qy_fractional,
        tune_knob1_name, tune_knob2_name,
        sequence_name, skip_use):
    print('\n Start coupling measurement...')
    if not skip_use:
        mad.use(sequence_name)

    # Store present values of tune knobs
    init_value = {}
    for kk in [tune_knob1_name, tune_knob2_name]:
        try:
            init_value[kk] = mad.globals[kk]
        except KeyError:
            print(f'{kk} not initialized, setting 0.0!')
            init_value[kk] = 0.


    # Try to push tunes on the diagonal
    qmid=(qx_fractional + qy_fractional)*0.5;
    qx_diagonal = qx_integer + qmid
    qy_diagonal = qy_integer + qmid
    mad.input(f'''
        match;
        global, q1={qx_diagonal},q2={qy_diagonal};
        vary,   name={tune_knob1_name}, step=1.E-9;
        vary,   name={tune_knob2_name}, step=1.E-9;
        lmdif,  calls=50, tolerance=1.E-10;
        endmatch;
    ''')

    # Measure closest tune approach
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    cta = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)


    # Restore intial values of tune knobs
    for kk in [tune_knob1_name, tune_knob2_name]:
        mad.globals[kk] = init_value[kk]

    print('\n Done coupling measurement.')

    return float(cta)


def coupling_correction(mad, n_iterations,
        qx_integer, qy_integer,
        qx_fractional, qy_fractional,
        tune_knob1_name, tune_knob2_name,
        cmr_knob_name, cmi_knob_name,
        sequence_name, skip_use):

    info_dict = {}
    print('\n Start coupling correction...')

    if not skip_use:
        mad.use(sequence_name)

    # Store present values of tune knobs
    init_value = {}
    for kk in [tune_knob1_name, tune_knob2_name]:
        try:
            init_value[kk] = mad.globals[kk]
        except KeyError:
            print(f'{kk} not initialized, setting 0.0!')
            init_value[kk] = 0.


    # Try to push tunes on the diagonal
    qmid=(qx_fractional + qy_fractional)*0.5;
    qx_diagonal = qx_integer + qmid
    qy_diagonal = qy_integer + qmid
    mad.input(f'''
        match;
        global, q1={qx_diagonal},q2={qy_diagonal};
        vary,   name={tune_knob1_name}, step=1.E-9;
        vary,   name={tune_knob2_name}, step=1.E-9;
        lmdif,  calls=50, tolerance=1.E-10;
        endmatch;
    ''')

    # Measure closest tune approach
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    cta0 = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)
    info_dict['closest_tune_appr_before_correction'] = cta0

    # Quick minimization based on linear machine
    cmrskew0 = mad.globals[cmr_knob_name]
    cmiskew0 = mad.globals[cmi_knob_name]
    info_dict['cmrknob_before_correction'] = cmrskew0
    info_dict['cmiknob_before_correction'] = cmiskew0

    # Optimization wrt to cmr
    mad.globals[cmr_knob_name] = cmrskew0 + cta0/2.;
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    ctap     = float(np.abs(2*(qx_tw-qy_tw) - np.round(2*(qx_tw-qy_tw)))/2)

    mad.globals[cmr_knob_name] = cmrskew0 - cta0/2.
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    ctam     = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)

    mad.globals[cmr_knob_name]= cmrskew0+(ctam**2-ctap**2)/2./cta0

    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    cta1     = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)


    # Optimization wrt to cmi
    mad.globals[cmi_knob_name] = cmiskew0 + cta1/2.
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    ctap     = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)

    mad.globals[cmi_knob_name]= cmiskew0-cta1/2.;
    mad.twiss()
    qx_tw= mad.table.summ.q1
    qy_tw= mad.table.summ.q2
    ctam     = float(np.abs(2*(qx_tw-qy_tw)-np.round(2*(qx_tw-qy_tw)))/2)

    mad.globals[cmi_knob_name]  = cmiskew0+(ctam**2-ctap**2)/2./cta1;


    #Empirical minimisation
    for i_iter in range(n_iterations):
        # Tolerances changed as in legacy routine
        if i_iter == 0:
            exp_tol = 6
        else:
            exp_tol = 7

        mad.input(f'''
            match;
            global, q1={qx_diagonal}, q2={qy_diagonal};
            vary,   name={tune_knob1_name}, step=1.E-9;
            vary,   name={tune_knob2_name}, step=1.E-9;
            lmdif,  calls=200,tolerance=1.E-{exp_tol};
            endmatch;

            match;
            global, q1={qx_diagonal}, q2={qy_diagonal};
            vary,   name={cmr_knob_name}, step=1.E-9;
            vary,   name={cmi_knob_name}, step=1.E-9;
            simplex,  calls=300, tolerance=2.E-{exp_tol};
            endmatch;
        ''')


    cta_after = coupling_measurement(mad, qx_integer, qy_integer,
        qx_fractional, qy_fractional,
        tune_knob1_name, tune_knob2_name,
        sequence_name, skip_use)
    info_dict['closest_tune_appr_after_correction'] = cta_after
    info_dict['cmrknob_after_correction'] = mad.globals[cmr_knob_name]
    info_dict['cmiknob_after_correction'] = mad.globals[cmi_knob_name]

    # Restore intial values of tune knobs
    for kk in [tune_knob1_name, tune_knob2_name]:
        mad.globals[kk] = init_value[kk]

    print('\n Done coupling correction.')

    return info_dict

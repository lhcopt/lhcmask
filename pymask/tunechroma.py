
def match_tune(mad, q1, q2,
        tune_knob1_name, tune_knob2_name,
        sequence_name, skip_use):

    if not skip_use:
        mad.use(sequence_name)

    mad.input(f'''
        match;
        global, q1={q1}, q2={q2};
        vary,   name={tune_knob1_name}, step=1.0E-7 ;
        vary,   name={tune_knob2_name}, step=1.0E-7 ;
        lmdif,  calls=100, tolerance=1.0E-21;
        endmatch;
        ''')

def match_chromaticity(mad, dq1, dq2,
        chromaticity_knob1_name, chromaticity_knob2_name,
        sequence_name, skip_use):

    if not skip_use:
        mad.use(sequence_name)

    mad.input(f'''
        match, chrom;
        global, dq1={dq1}, dq2={dq2};
        vary,   name={chromaticity_knob1_name};
        vary,   name={chromaticity_knob2_name};
        lmdif,  calls=100, tolerance=1.0E-21;
        endmatch;
        ''')

def match_tune_and_chromaticity(mad, q1, q2, dq1, dq2,
        tune_knob1_name, tune_knob2_name,
        chromaticity_knob1_name, chromaticity_knob2_name,
        sequence_name, skip_use):

    if not skip_use:
        mad.use(sequence_name)

    match_tune(mad, q1, q2,
        tune_knob1_name, tune_knob2_name,
        sequence_name, skip_use)

    match_chromaticity(mad, dq1, dq2,
        chromaticity_knob1_name, chromaticity_knob2_name,
        sequence_name, skip_use)

    mad.input(f'''
        match,chrom;
        global, dq1={dq1}, dq2={dq2};
        global, q1={q1}, q2={q2};
        vary,   name={chromaticity_knob1_name};
        vary,   name={chromaticity_knob2_name};
        vary,   name={tune_knob1_name}, step=1.0E-7 ;
        vary,   name={tune_knob2_name}, step=1.0E-7 ;
        lmdif,  calls=500, tolerance=1.0E-21;
        endmatch;
        ''')

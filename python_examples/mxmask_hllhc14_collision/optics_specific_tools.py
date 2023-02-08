import numpy as np

import pymask as pm

# The parts marked by (*) in following need to be
# adapted according to knob definitions

def build_sequence(mad, beam, **kwargs):

    assert 'optics_version' in kwargs.keys(), 'optics_version not specified'
    optics_version = kwargs['optics_version']

    # Select beam
    mad.input(f'mylhcbeam = {beam}')

    # Make link to optics toolkit
    pm.make_links({'optics_toolkit':
                f'optics_repository/toolkit'},
                force=True)

    mad.input('''

        ! Specify machine version
        ver_lhc_run = 0;
        '''
        f'''ver_hllhc_optics = {optics_version};'''
        f'''
        ! Get the toolkit
        call,file=
          "optics_repository/toolkit/macro.madx";
        '''
        '''
        ! Build sequence
        option, -echo,-warn,-info;
        if (mylhcbeam==4){
          call,file="optics_repository/../runIII/lhcb4.seq";
        } else {
          call,file="optics_repository/../runIII/lhc.seq";
        };
        option, -echo, warn,-info;
        '''
        f'''
        !Install HL-LHC
        call, file=
          "optics_repository/hllhc_sequence.madx";
        '''
        '''
        ! Slice nominal sequence
        exec, myslice;

        ! Install placeholder elements for errors (set to zero)
        call, file="errors/HL-LHC/install_MQXF_fringenl.madx";    ! adding fringe place holder
        call, file="errors/HL-LHC/install_MCBXFAB_errors.madx";   ! adding D1 corrector placeholders in IR1/5 (for errors)
        call, file="errors/HL-LHC/install_MCBRD_errors.madx";     ! adding D2 corrector placeholders in IR1/5 (for errors)
        call, file="errors/HL-LHC/install_NLC_errors.madx";       ! adding non-linear corrector placeholders in IR1/5 (for errors)

        !Cycling w.r.t. to IP3 (mandatory to find closed orbit in collision in the presence of errors)
        if (mylhcbeam<3){
          seqedit, sequence=lhcb1; flatten; cycle, start=IP3; flatten; endedit;
        };
        seqedit, sequence=lhcb2; flatten; cycle, start=IP3; flatten; endedit;

        ! Install crab cavities (they are off)
        call, file='optics_toolkit/enable_crabcavities.madx';
        on_crab1 = 0;
        on_crab5 = 0;

        ! Set twiss formats for MAD-X parts (macro from opt. toolkit)
        exec, twiss_opt;


        ''')

def apply_optics(mad, optics_file):
    mad.call(optics_file)


def set_optics_specific_knobs(mad, knob_settings, mode=None):

    # Copy knob settings to mad variable space
    mad.set_variables_from_dict(params=knob_settings)

    # A check
    if mad.globals.nrj < 500:
        assert knob_settings['on_disp'] == 0

    # A knob redefinition
    mad.input('on_alice := on_alice_normalized * 7000./nrj;')
    mad.input('on_lhcb := on_lhcb_normalized * 7000./nrj;')


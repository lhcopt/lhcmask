python_parameters = {

    # Links to be made for tools and scripts
    'links'                    : {
                                    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
                                    'modules': 'tracking_tools/modules',
                                    'tools': 'tracking_tools/tools',
                                    'beambeam_macros': 'tracking_tools/beambeam_macros',
                                    'errors': 'tracking_tools/errors',
                                 },
    # Mode - choose between:

    #   Main modes:
    #    'b1_without_bb'
    #    'b1_with_bb'
    #    'b4_from_b2_without_bb'
    #    'b4_from_b2_with_bb'

    #   Legacy modes
    #    'b1_with_bb_legacy_macros'
    #    'b4_without_bb'

    'mode'                      : 'b1_with_bb_legacy_macros',

    # Force separation in IP2 and IP8 if needed
    'force_leveling'            : None,

    # For test against madx mask for modes 'b4_from_b2_without_bb' and 'b4_without_bb':
    #'force_leveling' : {'on_sep8': -0.03425547139366354, 'on_sep2': 0.14471680504084292},

    # Tolerances for checks [ip1, ip2, ip5, ip8]
    'tol_beta'                 : [1e-3, 10e-2, 1e-3, 1e-2],
    'tol_sep'                  : [1e-6, 1e-6, 1e-6, 1e-6],

    # Tolerance for check on flat machine
    'tol_co_flatness'          : 1e-6,


    # Optics file
    'optics_file'              : '/afs/cern.ch/eng/lhc/optics/HLLHCV1.4/round/opt_round_150_1500_thin.madx', #15 cm

    # Enable checks
    'check_betas_at_ips'       : True,
    'check_separations_at_ips' : True,
    'save_intermediate_twiss'  : True,

    # Luminosity control in IP2 and IP8
    'enable_lumi_control'      : True,
    'sep_plane_ip2'            : 'x', # used by python tools - NOT by legacy macros
    'sep_plane_ip8'            : 'y', # used by python tools - NOT by legacy macros


    # Beam-beam parameters (used by python tools - NOT by legacy macros)
    'numberOfLRPerIRSide'      : [25, 20, 25, 20],
    'bunch_spacing_buckets'    : 10,
    'numberOfHOSlices'         : 11,
    'bunch_population_ppb'     : None,
    'sigmaz_m'                 : None,
    'z_crab_twiss'             : 0.075,

    # Enable multipolar errors
    'enable_multipolar_errors' : True,

    # N. iterations coupling correction
    'N_iter_coupling'            : 2,
    }


mask_parameters = {
    'par_verbose'              : 1,

    # Beam parameters
    'par_beam_norm_emit'       : 2.5,   # [um]
    'par_beam_sigt'            : 0.076,        # [m]
    'par_beam_sige'            : 1.1e-4,       # [-]
    'par_beam_npart'           : 2.2e11,       # [-]
    'par_beam_energy_tot'      : 7000,         # [GeV]

    # Settings
    'par_oct_current'          : -235,         # [A]
    'par_chromaticity'         : 5,            # [-] (Q':5 for colliding bunches, Q':15 for non-colliding bunches)
    'par_vrf_total'            : 16.,          # [MV]

    # Tunes
    'par_qx0'                  : 62.31,
    'par_qy0'                  : 60.32,


    #*************************#
    # Beam-beam configuration #
    #*************************#

    'par_on_bb_switch'         : 1,
    'par_match_with_bb'        : 0,            # should be off at collision
    'par_b_t_dist'             : 25.,          # bunch separation [ns]
    'par_n_inside_D1'          : 5,            # default value for the number of additionnal parasitic encounters inside D1

    'par_nho_IR1'              : 11,           # number of slices for head-on in IR1 (between 0 and 201)
    'par_nho_IR2'              : 11,           # number of slices for head-on in IR2 (between 0 and 201)
    'par_nho_IR5'              : 11,           # number of slices for head-on in IR5 (between 0 and 201)
    'par_nho_IR8'              : 11,           # number of slices for head-on in IR8 (between 0 and 201)

    #*****************************#
    #     Luminosity parameters   #
    #*****************************#

    # This variables set the leveled luminosity in IP8 
    'par_lumi_ip8'             : 2e33,         #[Hz/cm2]

    # This variables set the leveled luminosity in IP8 
    'par_fullsep_in_sigmas_ip2': 5,

    # These variables define the number of Head-On collisions in the 4 IPs
    'par_nco_IP1'              : 2748,
    'par_nco_IP2'              : 2494,
    'par_nco_IP5'              : 2748,
    'par_nco_IP8'              : 2572,

    #*************************#
    #  Errors and corrections #
    #*************************#

    # Select seed for errors
    'par_myseed'               : 1,

    # Set this flag to correct the errors of D2 in the NLC (warning: for now only correcting b3 of D2, still in development)
    'par_correct_for_D2'       : 0,
    # Set this flag to correct the errors of MCBXF in the NLC (warning: this might be less reproducable in reality, use with care)
    'par_correct_for_MCBX'     : 0,

    'par_on_errors_LHC'        : 1,
    'par_on_errors_MBH'        : 1,
    'par_on_errors_Q5'         : 1,
    'par_on_errors_Q4'         : 1,
    'par_on_errors_D2'         : 1,
    'par_on_errors_D1'         : 1,
    'par_on_errors_IT'         : 1,
    'par_on_errors_MCBRD'      : 0,
    'par_on_errors_MCBXF'      : 0,

}

knob_names = {
        # Common knobs
        'sepknob_ip2_mm': 'on_sep2',
        'sepknob_ip8_mm': 'on_sep8',

        # Knobs associated to sequences
        'qknob_1': {'lhcb1': 'kqtf.b1',  'lhcb2':'kqtf.b2'},
        'qknob_2': {'lhcb1': 'kqtd.b1',  'lhcb2':'kqtd.b2'},
        'chromknob_1': {'lhcb1': 'ksf.b1',  'lhcb2':'ksf.b2'},
        'chromknob_2': {'lhcb1': 'ksd.b1',  'lhcb2':'ksd.b2'},
        'cmrknob': {'lhcb1': 'cmrskew',  'lhcb2':'cmrskew'},
        'cmiknob': {'lhcb1': 'cmiskew',  'lhcb2':'cmiskew'},
        }

knob_parameters = {
    #IP specific orbit settings

    'par_x1'                   : 250,        # [urad]  # Could be optics-dependent
    'par_sep1'                 : 0,            # [mm]  # Could be optics-dependent
    'par_x2'                   : -170,         # [urad]  # Could be optics-dependent
    'par_sep2'                 : 0.138,        # [mm]  # Could be optics-dependent
    'par_x5'                   : 250,       # [urad]  # Could be optics-dependent
    'par_sep5'                 : 0,            # [mm]  # Could be optics-dependent
    'par_x8'                   : -250,         # [urad]  # Could be optics-dependent
    'par_sep8'                 : -0.043,       # [mm]  # Could be optics-dependent
    'par_a1'                   : 0,            # [urad]  # Could be optics-dependent
    'par_o1'                   : 0,            # [mm]  # Could be optics-dependent
    'par_a2'                   : 0,            # [urad]  # Could be optics-dependent
    'par_o2'                   : 0,            # [mm]  # Could be optics-dependent
    'par_a5'                   : 0,            # [urad]  # Could be optics-dependent
    'par_o5'                   : 0,            # [mm]  # Could be optics-dependent
    'par_a8'                   : 0,            # [urad]  # Could be optics-dependent
    'par_o8'                   : 0,            # [mm]  # Could be optics-dependent
    'par_crab1'                : -190,         # [urad]  # Could be optics-dependent
    'par_crab5'                : -190,    # [urad]  # Could be optics-dependent

    # Dispersion correction knob
    'par_on_disp'              : 1,            # Could be optics-dependent

    # Magnets of the experiments
    'par_on_alice'             : 1,
    'par_on_lhcb'              : 1,

    'par_on_sol_atlas'         : 0,
    'par_on_sol_cms'           : 0,
    'par_on_sol_alice'         : 0,
}

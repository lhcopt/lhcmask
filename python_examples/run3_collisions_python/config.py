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

    'mode'                      : 'b1_with_bb',

    # Force separation in IP2 and IP8 if needed
    'force_leveling'            : None,

    # Tolerances for checks [ip1, ip2, ip5, ip8]
    'tol_beta'                 : [1e-3, 10e-2, 1e-3, 1e-2],
    'tol_sep'                  : [1e-6, 1e-6, 1e-6, 1e-6],

    # Tolerance for check on flat machine
    'tol_co_flatness'          : 1e-6,

    # Optics file
    'optics_file'              : '/afs/cern.ch/eng/lhc/optics/runIII/RunIII_dev/2022_V1/PROTON/opticsfile.29',

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
    'z_crab_twiss'             : 0.,

    # Enable machine imperfections
    'enable_imperfections'     : False,

    # Enable knob synthesis (for coupling correction, if no errors)
    'enable_knob_synthesis'    : False,

    # Enable crab cavities
    'enable_crabs'             : False,

    # N. iterations coupling correction
    'N_iter_coupling'            : 2,

    # Value to be added to linear coupling knobs (on sequence_to_track)
    'delta_cmr'                 : 1e-3,
    'delta_cmi'                 : 0.,
    }


mask_parameters = {
    'par_verbose'              : 1,

    # Beam parameters
    'par_beam_norm_emit_x'     : 2.5,          # [um]
    'par_beam_norm_emit_y'     : 2.5,          # [um]
    'par_beam_sigt'            : 0.076,        # [m]
    'par_beam_sige'            : 1.1e-4,       # [-]
    'par_beam_npart'           : 1.8e11,       # [-]
    'par_beam_energy_tot'      : 7000,         # [GeV]

    # Settings
    'par_oct_current'          : -350,         # [A]
    'par_chromaticity_x'       : 15,            # [-] 
    'par_chromaticity_y'       : 15,           # [-] 
    'par_vrf_total'            : 12.,          # [MV]

    # Tunes
    'par_qx0'                  : 62.313,
    'par_qy0'                  : 60.318,


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
    'par_nco_IP1'              : 2736,
    'par_nco_IP2'              : 2250,
    'par_nco_IP5'              : 2736,
    'par_nco_IP8'              : 2376,

    #*************************#
    #  Errors and corrections #
    #*************************#

    # Select seed for errors
    'par_myseed'               : 0,

    # Set this flag to correct the errors of D2 in the NLC (warning: for now only correcting b3 of D2, still in development)
    'par_correct_for_D2'       : 0,
    # Set this flag to correct the errors of MCBXF in the NLC (warning: this might be less reproducable in reality, use with care)
    'par_correct_for_MCBX'     : 0,

    'par_on_errors_LHC'        : 0,
    'par_on_errors_MBH'        : 0,
    'par_on_errors_Q5'         : 0,
    'par_on_errors_Q4'         : 0,
    'par_on_errors_D2'         : 0,
    'par_on_errors_D1'         : 0,
    'par_on_errors_IT'         : 0,
    'par_on_errors_MCBRD'      : 0,
    'par_on_errors_MCBXF'      : 0,

}

knob_names = {
        # Common knobs
        'sepknob_ip2_mm': 'on_sep2h',
        'sepknob_ip8_mm': 'on_sep8v',

        # Knobs associated to sequences
        'qknob_1': {'lhcb1': 'dQx.b1_sq',  'lhcb2':'dQx.b2_sq'},
        'qknob_2': {'lhcb1': 'dQy.b1_sq',  'lhcb2':'dQy.b2_sq'},
        'chromknob_1': {'lhcb1': 'dQpx.b1_sq',  'lhcb2':'dQpx.b2_sq'},
        'chromknob_2': {'lhcb1': 'dQpy.b1_sq',  'lhcb2':'dQpy.b2_sq'},
        'cmrknob': {'lhcb1': 'CMRS.b1_sq',  'lhcb2':'CMRS.b2_sq'},
        'cmiknob': {'lhcb1': 'CMIS.b1_sq',  'lhcb2':'CMIS.b2_sq'},
        }

knob_settings = {
    #IP specific orbit settings

    'on_x1'                   : 150,          # [urad]  
    'on_sep1'                 : 0,            # [mm]   
    'on_x2h'                  : 0,            # [urad] 
    'on_x2v'                  : 200,          # [urad] 
    'on_sep2h'                : 1.0,          # [mm]   
    'on_sep2v'                : 0,            # [mm]   
    'on_x5'                   : 150,          # [urad] 
    'on_sep5'                 : 0,            # [mm]   
    'on_x8h'                  : -250,         # [urad] 
    'on_x8v'                  : 0,            # [urad] 
    'on_sep8h'                : 0.,           # [mm]   
    'on_sep8v'                : -1.0,         # [mm]   

    # Dispersion correction knob
    'on_disp'                 : 0.,           # Could be optics-dependent

    # Magnets of the experiments
    'on_alice_normalized'     : 1,
    'on_lhcb_normalized'      : 1,

    'on_sol_atlas'            : 0,
    'on_sol_cms'              : 0,
    'on_sol_alice'            : 0,
}

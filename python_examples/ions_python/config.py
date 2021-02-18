configuration = {

    # Links to be made for tools and scripts
    'links'                    : {
                                    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
                                    #'modules': 'tracking_tools/modules',
                                    'modules': '../../',
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

    # Optics file
    'optics_file'              : '/afs/cern.ch/eng/lhc/optics/runII/2018/ION/opticsfile.21',

    # Enable checks
    'check_betas_at_ips'       : True,
    'check_separations_at_ips' : True,
    'save_intermediate_twiss'  : True,

   # Tolerances for checks [ip1, ip2, ip5, ip8]
    'tol_beta'                 : [10e-2, 10e-2, 10e-2, 20e-2],
    'tol_sep'                  : [1e-6, 1e-6, 1e-6, 1e-6],

    # Tolerance for check on flat machine
    'tol_co_flatness'          : 1e-6,

    # Beam parameters
    'beam_norm_emit_x'     : 1.65,          # [um]
    'beam_norm_emit_y'     : 1.65,          # [um]
    'beam_sigt'            : 0.0824,        # [m]
    'beam_sige'            : 1.02e-4,       # [-]
    'beam_npart'           : 1.8e8,         # [-]
    'beam_energy_tot'      : 574000.,       # [GeV]

    # Ion parameters
    'mass'                 :193.68715,
    'charge'               :82,

    # Tunes and chromaticities
    'qx0'                  : 62.31,
    'qy0'                  : 60.32,
    'chromaticity_x'       : 10,            # [-] 
    'chromaticity_y'       : 10,            # [-] 

    # RF voltage
    'vrf_total'            : 14.,          # [MV]

    # Octupole current
    'oct_current'          : 250,         # [A]

    # Luminosity parameters
    'enable_lumi_control'      : True,
    'sep_plane_ip1'            : 'x', # used by python tools - NOT by legacy macros
    'sep_plane_ip2'            : 'x', # used by python tools - NOT by legacy macros
    'sep_plane_ip5'            : 'y', # used by python tools - NOT by legacy macros
    'sep_plane_ip8'            : 'y', # used by python tools - NOT by legacy macros
    'lumi_ip8'             : 1e27,           #[Hz/cm2]
    'lumi_ip1'             : 6.4e27,         #[Hz/cm2]
    'lumi_ip2'             : 6.4e27,         #[Hz/cm2]
    'lumi_ip5'             : 6.4e27,         #[Hz/cm2]
    #'fullsep_in_sigmas_ip2': 5,
    'nco_IP1'              : 1088,
    'nco_IP2'              : 1088,
    'nco_IP5'              : 1088,
    'nco_IP8'              : 398,

    # Beam-beam parameters (used by python tools - NOT by legacy macros)
    'numberOfLRPerIRSide'      : [25, 20, 25, 20],
    'bunch_spacing_buckets'    : 20,
    'numberOfHOSlices'         : 11,
    'bunch_population_ppb'     : None,
    'sigmaz_m'                 : None,
    'z_crab_twiss'             : 0.,

    # Match tunes and chromaticities including beam-beam effects
    'match_q_dq_with_bb'        : False,            # should be off at collision

    # Enable crab cavities
    'enable_crabs'             : False,

    # N. iterations coupling correction
    'N_iter_coupling'            : 2,

    # Value to be added to linear coupling knobs (on sequence_to_track)
    'delta_cmr'                 : 0.,
    'delta_cmi'                 : 0.,

    # Verbose flag for MAD-X parts
    'verbose_mad_parts'         : True,

    # Optics-specific knob namings
    'knob_names' : {
        # Common knobs
        'sepknob_ip2_mm': 'on_sep2',
        'sepknob_ip8_mm': 'on_sep8',
        'sepknob_ip1_mm': 'on_sep1',
        'sepknob_ip5_mm': 'on_sep5',

        # Knobs associated to sequences
        'qknob_1': {'lhcb1': 'dQx.b1_sq',  'lhcb2':'dQx.b2_sq'},
        'qknob_2': {'lhcb1': 'dQy.b1_sq',  'lhcb2':'dQy.b2_sq'},
        'chromknob_1': {'lhcb1': 'dQpx.b1_sq',  'lhcb2':'dQpx.b2_sq'},
        'chromknob_2': {'lhcb1': 'dQpy.b1_sq',  'lhcb2':'dQpy.b2_sq'},
        'cmrknob': {'lhcb1': 'CMRS.b1_sq',  'lhcb2':'CMRS.b2_sq'},
        'cmiknob': {'lhcb1': 'CMIS.b1_sq',  'lhcb2':'CMIS.b2_sq'},
        },

    # Optics specific knob values
    # (only on_disp is used directly by the mask,
    # the other values are used only throught the optics_specific_tools file)
    'knob_settings':  {
        #IP specific orbit settings
        'on_x1'                   : 170,          # [urad]  
        'on_sep1'                 : 1e-8,         # [mm]   
        'on_x2'                   : 170,          # [urad] 
        'on_sep2'                 : 1e-8,         # [mm]   
        'on_x5'                   : 170,          # [urad] 
        'on_sep5'                 : 1e-8,         # [mm]   
        'on_x8'                   : -170,         # [urad] 
        'on_sep8'                 : 1e-8,         # [mm]   
        'on_ov2'                  :0.,
        'on_ov5'                  :0.,
        # Dispersion correction knob
        'on_disp'                 : 0,            # Value to choose could be optics-dependent

        # Magnets of the experiments
        'on_alice_normalized'     : -1,
        'on_lhcb_normalized'      : -1,

        'on_sol_atlas'            : 0,
        'on_sol_cms'              : 0,
        'on_sol_alice'            : 0,
        },

    # Enable machine imperfections
    'enable_imperfections'     : False,

    # Enable knob synthesis (for coupling correction, if no imperfections)
    'enable_knob_synthesis'    : False,

    # Parameters for machine imperfections
    'pars_for_imperfections': {},

    # Parameters for legacy beam-beam macros (not used in default modes)
    'pars_for_legacy_bb_macros' : {
                                    'par_b_t_dist' : 50.,  # bunch spacing [ns]
                                    'par_n_inside_D1': 5,  # n. parasitic encounters inside D1
                                  },
    }

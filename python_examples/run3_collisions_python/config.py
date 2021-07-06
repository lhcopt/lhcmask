configuration = {

    # Links to be made for tools and scripts
    'links'                    : {
                                    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
                                    'modules': 'tracking_tools/modules',
                                    'tools': 'tracking_tools/tools',
                                    'beambeam_macros': 'tracking_tools/beambeam_macros',
                                    'errors': 'tracking_tools/errors',
                                    'optics_repository' : '/afs/cern.ch/eng/lhc/optics',
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
    'optics_file'              : 'optics_repository/runIII/RunIII_dev/2022_V1/PROTON/opticsfile.29',

    # Enable checks
    'check_betas_at_ips'       : True,
    'check_separations_at_ips' : True,
    'save_intermediate_twiss'  : True,

   # Tolerances for checks [ip1, ip2, ip5, ip8]
    'tol_beta'                 : [1e-3, 10e-2, 1e-3, 1e-2],
    'tol_sep'                  : [1e-6, 1e-6, 1e-6, 1e-6],

    # Tolerance for check on flat machine
    'tol_co_flatness'          : 1e-6,

    # Beam parameters
    'beam_norm_emit_x'     : 2.5,          # [um]
    'beam_norm_emit_y'     : 2.5,          # [um]
    'beam_sigt'            : 0.076,        # [m]
    'beam_sige'            : 1.1e-4,       # [-]
    'beam_npart'           : 1.8e11,       # [-]
    'beam_energy_tot'      : 7000,         # [GeV]

    # Tunes and chromaticities
    'qx0'                  : 62.313,
    'qy0'                  : 60.318,
    'chromaticity_x'       : 15,            # [-] 
    'chromaticity_y'       : 15,            # [-] 

    # RF voltage
    'vrf_total'            : 12.,          # [MV]

    # Octupole current
    'oct_current'          : -350,         # [A]

    # Luminosity parameters
    'enable_lumi_control'      : True,
    'sep_plane_ip2'            : 'x', # used by python tools - NOT by legacy macros
    'sep_plane_ip8'            : 'y', # used by python tools - NOT by legacy macros
    'lumi_ip8'             : 2e33,         #[Hz/cm2]
    'fullsep_in_sigmas_ip2': 5,
    'nco_IP1'              : 2736,
    'nco_IP2'              : 2250,
    'nco_IP5'              : 2736,
    'nco_IP8'              : 2376,

    # Beam-beam parameters (used by python tools - NOT by legacy macros)
    'beambeam_config'      :
        {
            'numberOfLRPerIRSide'      : [25, 20, 25, 20],
            'bunch_spacing_buckets'    : 10,
            'numberOfHOSlices'         : 11,
            'bunch_num_particles'      : None,
            'bunch_particle_charge'    : None,
            'sigmaz_m'                 : None,
            'z_crab_twiss'             : 0.,
        },

    # Match tunes and chromaticities including beam-beam effects
    'match_q_dq_with_bb'        : False,            # should be off at collision

    # Enable crab cavities
    'enable_crabs'             : False,

    # N. iterations for coupling correction
    'N_iter_coupling'            : 2,

    # Value to be added to linear coupling knobs (on sequence_to_track)
    'delta_cmr'                 : 1e-3,
    'delta_cmi'                 : 0.,

    # Verbose flag for MAD-X parts
    'verbose_mad_parts'         : True,

    # Optics-specific knob namings
    'knob_names' : {
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
        },

    # Optics specific knob values
    # (only on_disp is used directly by the mask,
    # the other values are used only throught the optics_specific_tools file)
    'knob_settings':  {
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
        'on_disp'                 : 0,            # Value to choose could be optics-dependent

        # Magnets of the experiments
        'on_alice_normalized'     : 1,
        'on_lhcb_normalized'      : 1,

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
                                    'par_b_t_dist' : 25.,  # bunch spacing [ns]
                                    'par_n_inside_D1': 5,  # n. parasitic encounters inside D1
                                  },
    }

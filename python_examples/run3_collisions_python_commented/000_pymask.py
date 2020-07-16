# %%
'''
# A pythonic approach to run LHC mask for Run 3

In this script, we are going to describe, step-by-step, the 
approach we propose for  Run 3 masks (and, more in general,for all LHC masks) in python.

The main purpose of the new approach is to profit 

1. from the scripting capability of python and
2. from beam optics computation of MAD-X.

!!! info
	Our aim is to move as much as possible of the "scripting" logic from MAD-X
	to python, as we consider the latter much more flexible for scripting purpose.
	This "translation" is still on-going and will require sometime. Hereby we focus
	on define the correct working-flow based by defining interfaces between MAD-X 
	and python..

The first application we focus on is the beam-beam, BB, studies of B4 
(not possible with the legacy approach). 

To concentrate on the contents and not to waste time to set up on the environment, 
you can login on a *lxplus* machine and do 

```bash
source /afs/cern.ch/eng/tracking-tools/python_installations/activate_default_python  
```

This sets up a python installation with the required packages and tools.

Then you can clone the *lhcmask* repository on your terminal with

```bash
cd
git clone https://github.com/lhcopt/lhcmask.git
```

then you just go on the folder 

```bash
cd ~/lhcmask/python_examples/run3_collisions_python_commented
```

and the file we are commenting is the **000_pymask.py**.

!!! warning
	The previous folder does not exist yet in the master repository. 
	For the moment you have only the folder *run3_collisions_python* 
	(i.e., the version without comments).

### Usual import of packages 
'''

# %%
import sys, os, pickle, numpy as np

# %%
'''
As we will see, in this script we will use 

- optics dependent and
- optics independent

modules/packages.

The *pymask* package is optics independent and comes with the *lhcmask*
repository. It contains

- the Madxp utility to access the MAD-X constant or dependent/independent variables 
- *make_links* to link the different needed folders on the local directory.
- *checks_on_parameter_dict* 
- *get_pymask_configuration*
- *generate_bb_dataframes*
- *configure_b4_from_b2*
- *install_lenses_in_sequence*
- *generate_sixtrack_input*
- *get_optics_and_orbit_at_start_ring*
- *generate_pysixtrack_line_with_bb*
- ...

!!! info
	The user defines and takes care of the optics dependent 
	modules/packages. He can edit the content of this file (the mask file). 
	On the other hand, all the MAD-X modules should not be edited by the user and 
	the user should point to the repository
	
	```bash
	/afs/cern.ch/eng/tracking-tools/...
	```
  
	This will avoid code duplication and to maintain the code in a centralized 
  repository.
	All folders in **tracking-tools** are git repositories
  (all but the **python_installations**).
	The user (to control and freeze at a given point the repository) can clone them
  locally BUT then is her/his responsibility to check systematically
  if the local git repositories are up-to-date with the master repository.
'''

# %%
sys.path.append('/afs/cern.ch/eng/tracking-tools/modules/')
import pymask as pm
Madx = pm.Madxp

# %%
'''The user has the possibility to provide some optics specific functions and 
dictionaries of parameters and knobs. Namely:

- optics_specific_tools.py: a set of custom functions.
- mask_parameters.py: the dictionary of the parameters to be used in 
the simulations. 
- knob_parameters.py: the dictionary of the knobs that set the sequence(s).

!!! info
	A knob of a sequence is a independent variables of MAD-X that changes
	the attributes of the element of the sequence(s). E.g., $on_x1$ is a knob 
	but $qx0$ is a parameter.

The three files are located in the working folder. 
In the optics specific tools we will find methods like:

- *build_sequence*
- *apply_optics*
- *twiss_and_check*
- ...
'''

# %%
import optics_specific_tools as ost
from mask_parameters import mask_parameters

# %%
'''
### Selecting the beam mode
The user needs then to select the beam mode. There are 
6 beam modes presently available:

- b1_without_bb
- b1_with_bb
- b1_with_bb_legacy_macros
- b4_without_bb
- b4_from_b2_without_bb
- b4_from_b2_with_bb

Their purpose is quite explicit in their names. 
It is important to note that the B4 has not legacy BB macros (since it is not possible
to run the B4 with the legacy macros).
'''

# %%
mode = 'b4_from_b2_with_bb'

# %%
'''
The flexibility in python allows us to make checks during the mask in a more 
systematic way. The users should introduce them for monitoring
the sanity of the code.

As we will show, there are optics independent function that help the assertions.
The users needs to define some tolerances, e.g., IP dependent tolerances for the
beta-function or the beam separations.
'''

# %%
tol_beta = [1e-3, 10e-2, 1e-3, 1e-2]
tol_sep = [1e-6, 1e-6, 1e-6, 1e-6]

# %%
'''
### Defining the links
To define the links to the repository folders one can use the **pymask** function.
As mentioned, all links should point to /afs/cern.ch/eng/tracking-tools.

!!! tip
	Again: it is important to maintaing a central version of the MAD-X modules and macros 
	to ease their maintenance and to ensure the natural propagation of the fixes. 
	So the users should customize this script and the optics-dependent modules but 
	not the optics independent ones. Clearly is more than welcome to contribure by 
	pulling a request in the git repository.
'''

# %%
pm.make_links(force=True, links_dict={
    'tracking_tools': '/afs/cern.ch/eng/tracking-tools',
    'modules': 'tracking_tools/modules',
    'tools': 'tracking_tools/tools',
    'beambeam_macros': 'tracking_tools/beambeam_macros',
    'errors': 'tracking_tools/errors'})

# %%
'''
Define few additional parameters.
'''

# %%
optics_file = 'opticsfile.29'
# to check the beta at the ips
check_betas_at_ips = True
# to check the separation at the ips
check_separations_at_ips = True
# to save in parquet files the intermediate twiss
# for further checks or plots 
save_intermediate_twiss = True

# %%
'''
### Check and load the parameters.
'''

# %%
pm.checks_on_parameter_dict(mask_parameters)

# %%
'''
### Configuration definition

Depending on the beam mode, different flags will be set accordingly.
In fact a beam mode is nothing else that a consistent set of flags.
'''

# %%
(beam_to_configure, sequences_to_check, sequence_to_track, generate_b4_from_b2,
    track_from_b4_mad_instance, enable_bb_python, enable_bb_legacy,
    force_disable_check_separations_at_ips,
    ) = pm.get_pymask_configuration(mode)

if force_disable_check_separations_at_ips:
    check_separations_at_ips = False

# %%
'''
### Starting MAD-X
'''

# %%
mad = Madx()

# %%
'''
###  Build the sequence

This is the first example of function that the user need to define and it is 
optics dependent.
'''

# %%
ost.build_sequence(mad, beam=beam_to_configure)

# %%
'''
### Load a specific optics

This is the second example of user's defined function.
'''

# %%
ost.apply_optics(mad, optics_file=optics_file)

# %%
'''
### Force disable beam-beam when needed
'''

# %%
if not(enable_bb_legacy) and not(enable_bb_python):
    mask_parameters['par_on_bb_switch'] = 0.

# %%
'''
### Load parameters to mad
'''

# %%
mad.set_variables_from_dict(params=mask_parameters)

# %%
'''
### Prepare sequences and attach beam
Now we start to call the MAD-X modules. They are optics independent.
'''

# %%
mad.call("modules/submodule_01a_preparation.madx")
mad.call("modules/submodule_01b_beam.madx")

# %%
'''
### SANITY CHECK: test the machine from the repository
Here is important to note that the user's knob are not yet applied. 
Therefore the machine has the knob set from the optics repository.
'''

# %%
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_from_optics',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips,
        check_separations_at_ips=check_separations_at_ips)

# %%
'''
### Set phase, apply and save crossing
'''

# %%
mad.call("modules/submodule_01c_phase.madx")

# %%
'''
### Set optics-specific knobs
Here the user needs to define the good conventions.
Some of the knobs definition are optics dependent (unfortunately). 
The user need to make the proper knob conversion.  
'''

# %%
ost.set_optics_specific_knobs(mad, mode)

# %%
'''
### Crossing-save and some reference measurements
'''

# %%
mad.input('exec, crossing_save')
mad.call("modules/submodule_01e_final.madx")

# %%
'''
### Test flat machine
'''

# %%
mad.input('exec, crossing_disable')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_no_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)

# %%
'''
### Check CO flatness
'''

# %%
flat_tol = 1e-6
for ss in twiss_dfs.keys():
    tt = twiss_dfs[ss]
    assert np.max(np.abs(tt.x)) < flat_tol
    assert np.max(np.abs(tt.y)) < flat_tol

# %%
'''
### Check machine after crossing restore
'''

# %%
mad.input('exec, crossing_restore')
twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_with_crossing',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=check_separations_at_ips)

mad.use(f'lhcb{beam_to_configure}')

# %%
'''
### Call leveling module
'''

# %%
# if mode=='b4_without_bb':
#     print('Leveling not working in this mode!')
# else:
#     mad.call("modules/module_02_lumilevel.madx")
mad.input('on_disp = 0')

# %%
'''
### Prepare bb dataframes
'''

# %%
if enable_bb_python:
    bb_dfs = pm.generate_bb_dataframes(mad,
        ip_names=['ip1', 'ip2', 'ip5', 'ip8'],
        harmonic_number=35640,
        numberOfLRPerIRSide=[25, 20, 25, 20],
        bunch_spacing_buckets=10,
        numberOfHOSlices=11,
        bunch_population_ppb=None,
        sigmaz_m=None,
        z_crab_twiss = 0,
        remove_dummy_lenses=True)

# %%
'''
!!! info
	Here the dataframes can be edited, e.g., to set bbb intensity
'''

# %%
'''
### Generate mad instance for B4
'''

# %%
if generate_b4_from_b2:
    mad_b4 = Madx()
    ost.build_sequence(mad_b4, beam=4)
    ost.apply_optics(mad_b4, optics_file=optics_file)

    pm.configure_b4_from_b2(mad_b4, mad)

    twiss_dfs_b2, other_data_b2 = ost.twiss_and_check(mad,
            sequences_to_check=['lhcb2'],
            tol_beta=tol_beta, tol_sep=tol_sep,
            twiss_fname='twiss_b2_for_b4check',
            save_twiss_files=save_intermediate_twiss,
            check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

    twiss_dfs_b4, other_data_b4 = ost.twiss_and_check(mad_b4,
            sequences_to_check=['lhcb2'],
            tol_beta=tol_beta, tol_sep=tol_sep,
            twiss_fname='twiss_b4_for_b4check',
            save_twiss_files=save_intermediate_twiss,
            check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

# %%
''' 
### From collider to synchrotron
We working exclusively on the sequence to track.
'''

# %%
# Select mad object
if track_from_b4_mad_instance:
    mad_track = mad_b4
else:
    mad_track = mad

mad_collider = mad
del(mad)

# %%
'''
### Twiss machine to track
'''

# %%
twiss_dfs, other_data = ost.twiss_and_check(mad_track, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_track_intermediate',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

# %%
'''
### Install the bb lenses
'''

# %%
if enable_bb_python:
    if track_from_b4_mad_instance:
        bb_df_track = bb_dfs['b4']
        assert(sequence_to_track=='lhcb2')
    else:
        bb_df_track = bb_dfs['b1']
        assert(sequence_to_track=='lhcb1')

    pm.install_lenses_in_sequence(mad_track, bb_df_track, sequence_to_track)

    # Disable bb
    mad_track.globals.on_bb_charge = 0
else:
    bb_df_track = None

# %%
'''
### Legacy bb macros
'''

# %%
if enable_bb_legacy:
    assert(beam_to_configure == 1)
    assert(not(track_from_b4_mad_instance))
    assert(not(enable_bb_python))
    mad_track.call("modules/module_03_beambeam.madx")

# %%
'''
### Install crab cavities
'''

# %%
# mad_track.call("tools/enable_crabcavities.madx")
# mad_track.call("modules/submodule_04_1a_install_crabs.madx")

# %%
'''
### Save references (orbit at IPs)
'''

# %%
mad_track.call('modules/submodule_04_1b_save_references.madx')

# %%
'''
### Switch off dipersion correction knob
'''

# %%
mad_track.globals.on_disp = 0.

# %%
'''
### Final use of "use"
'''
mad_track.use(sequence_to_track)

# %%
'''
#### Disabling "use"
'''

# %%
mad_track._use = mad_track.use
mad_track.use = None

# %%
'''
### Install and correct errors
'''

# %%
# mad_track.call("modules/module_04_errors.madx")

# %%
'''
### Machine tuning (enables bb)
'''

# %%
mad_track.call("modules/module_05_tuning.madx")

# %%
'''
### Switch on crab cavities
'''

# %%
# mad_track.globals.on_crab1 = mad_track.globals.par_crab1
# mad_track.globals.on_crab5 = mad_track.globals.par_crab5

# %%
'''
### Generate sixtrack input
'''

# %%
if enable_bb_legacy:
    mad_track.call("modules/module_06_generate.madx")
else:
    pm.generate_sixtrack_input(mad_track,
        seq_name=sequence_to_track,
        bb_df=bb_df_track,
        output_folder='./',
        reference_bunch_charge_sixtrack_ppb=(
            mad_track.sequence[sequence_to_track].beam.npart),
        emitnx_sixtrack_um=(
            mad_track.sequence[sequence_to_track].beam.exn),
        emitny_sixtrack_um=(
            mad_track.sequence[sequence_to_track].beam.eyn),
        sigz_sixtrack_m=(
            mad_track.sequence[sequence_to_track].beam.sigt),
        sige_sixtrack=(
            mad_track.sequence[sequence_to_track].beam.sige),
        ibeco_sixtrack=1,
        ibtyp_sixtrack=0,
        lhc_sixtrack=2,
        ibbc_sixtrack=0,
        radius_sixtrack_multip_conversion_mad=0.017,
        skip_mad_use=True)

# %%
'''
### Get optics and orbit at the start ring
'''

# %%
optics_orbit_start_ring = pm.get_optics_and_orbit_at_start_ring(
        mad_track, sequence_to_track, skip_mad_use=True)
with open('./optics_orbit_at_start_ring.pkl', 'wb') as fid:
    pickle.dump(optics_orbit_start_ring, fid)

# %%
'''
### Generate pysixtrack lines
'''

# %%
if enable_bb_legacy:
    print('Pysixtrack line is not generated with bb legacy macros')
else:
    pysix_fol_name = "./pysixtrack"
    dct_pysxt = pm.generate_pysixtrack_line_with_bb(mad_track,
        sequence_to_track, bb_df_track,
        closed_orbit_method='from_mad',
        pickle_lines_in_folder=pysix_fol_name,
        skip_mad_use=True)

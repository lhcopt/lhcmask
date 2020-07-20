# %%
'''
# A pythonic approach to run LHC mask for Run 3

In this script, we are going to describe, step-by-step, the 
approach we propose for Run 3 mask(s) in python.
This applies, more in general, to all LHC and HL-LHC masks.

!!! info
	Here we are using the legacy jargon of the "mask": in reality one could simply 
	refer to the "mask" as a (python) script. We will discuss here the structure of
	the script and the rationale we adopted.
	

The main purpose of the new approach is to profit 

1. from the scripting capability of python and
2. from beam optics computation of MAD-X.

For example, the luminosity leveling is nothing to do with the computing machine 
optics (MAD-X domain) and is much more natural to do it in python.

!!! info
	Our aim is to move as much as possible of the "scripting" logic from MAD-X
	to python, as we consider the latter much more flexible for scripting purpose.
	This "translation" is still on-going and will require sometime. Hereby we focus
	on define the correct working-flow based by defining interfaces between MAD-X 
	and python..

The first application we focus on is the beam-beam, BB, studies of B4.
This study is not possible with the legacy approach and is a good example
of the gain in the flexibility python can provide. 

To concentrate on the contents and not to spend time to setup the environment, 
you can login on a *lxplus* machine and execute 

```bash
source /afs/cern.ch/eng/tracking-tools/python_installations/activate_default_python  
```

This setups a python installation with the required packages and tools.

Then you can clone the *lhcmask* repository on your terminal with

```bash
cd
git clone https://github.com/sterbini/lhcmask.git
# I will make a pull request later so that you can do
# git clone https://github.com/lhcopt/lhcmask.git
# and use the official repository
```

then you just go on the folder 

```bash
cd ~/lhcmask/python_examples/run3_collisions_python_commented
```

and the file we are commenting is the **000_pymask.py**.

!!! warning
	The previous folder does not exist yet in the official repository. 
	For the moment you have it only in my fork.

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
repository (together with most of the optics independent modules). It contains

- the Madxp module to access the MAD-X variable workspace in a hierarchical way
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
	modules/packages. He can edit the content of this script (the mask file). 
	On the other hand, all the MAD-X modules should not be edited by the user and 
	the user should point to the repository
	
	```bash
	/afs/cern.ch/eng/tracking-tools/...
	```
  
	This will avoid code duplication and help to maintain the code in a centralized 
  repository.
	All folders in **tracking-tools** are git repositories
  (all but the **python_installations**, the one that contains the python
	distribution we suggested to source).
	
	Clearly we are aware that the user wants to have the full control of the code, 
	e.g., to ensure that there are no updates or new realeases in the middle of a simulation.
	
	To tackle the problem, the user can clone the corresponding git repositories locally
  BUT then is her/his responsibility to check systematically
  if her/his local git repositories are up-to-date with the master repository.

	The user is more than welcome to contribute by forking the repository and 
	pull-requesting her/his contribution.
'''

# %%
# sys.path.append('/afs/cern.ch/eng/tracking-tools/modules/')
sys.path.append('../../')
import pymask as pm
Madx = pm.Madxp
lumi = pm.luminosity

# %%
'''The user has to provide the optics specific functions and the 
dictionaries of parameters and knobs. Namely:

- optics_specific_tools.py: a set of custom functions.
- mask_parameters.py: the dictionary of the parameters to be used in 
the simulations. 
- knob_parameters.py: the dictionary of the knobs that set the sequence(s).

!!! info
	A knob of a sequence is a independent variables of MAD-X that changes
	the attributes of the element of the sequence. E.g., $on_x1$ is a knob 
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

A beam mode will set several flags. We operate in the "beam mode space" 
instead in the "flag space" since not all the flag combinations are physicals. 
'''

# %%
# here we consider the b4 built from the b2 and with BB.
mode = 'b4_from_b2_with_bb'

# %%
'''
!!! info
	We would like to make a digression on the B4. For periodic linear problem we can
	find the twiss the machine optics for B2 or B4 and find an equivalent solutions.
	The tracking is differnt: let's assume a kicker introducing a betatronic 
	oscillation: clearly the direction of the beam is important to describe the
	beam trajectory.
	
	In the optics repository we have B1, B2 and B4 but we miss B3. For the BB 
	we needs B1 and B2 or B3 and B4. Since we do not have B3, in order
	to have BB in B4, we compute the B1/B2 BB configuration. We apply the BB lenses
	in B2, then starting from the B2 elements we populates the attributes of the 
	B4 sequences (we call it the B2-to-B4 transplant).
'''

# %%
'''
The flexibility in python allows us to make checks during the mask in a more 
systematic way. The user should introduce them for monitoring
the sanity of the code.

As we will show, there are optics independent function that help the sanity 
assertions.
The user needs to define some tolerances, e.g., IP dependent tolerances for the
beta-function or the beam separations.
'''

# %%
tol_beta = [1e-3, 10e-2, 1e-3, 1e-2] # in m
tol_sep = [1e-6, 1e-6, 1e-6, 1e-6]	 # in m

# %%
'''
### Defining the links
To define the links to the repository folders one can use the **pymask** function.
As mentioned, all links should point to /afs/cern.ch/eng/tracking-tools 
(or to the clones of the relevant git repository).

!!! tip
	Again: we encourage to maintain a central version of the MAD-X modules and macros 
	to ease their maintenance and to ensure the natural propagation of the fixes. 
	So the users should customize this script and the optics-dependent modules but 
	not the optics independent ones. Clearly the user is more than welcome to
	contribute by forking, editing and pulling a request in the git repository.
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
Define few additional parameters used in the script. The user can custom the
script as she/he likes.
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

!!! tip
	Please use the python namespace to associate the fuction to a specific
	package/modules and use
	```python
	import inspect
	print(''.join(inspect.getsourcelines(pm.checks_on_parameter_dict)[0]))
	```
	to inspect the function code or to read the function help (and remember the
	difficulties to find similar information for MAD-X macros). On the same line,
	please use the python debugger (pdb) to debug your script.
	
	```bash
	python -m pdb 000_pymask.py 	
  ```
	
	and once you are in the pdb you can go to line 276 by
	
	```python
	break 276
	continue
  ```
	
	and then inspect the variables and debug the code.
'''

# %%
pm.checks_on_parameter_dict(mask_parameters)

# %%
'''
### Configuration definition

Depending on the beam mode, different flags will be set accordingly.
In fact, as already mentioned, a beam mode is nothing else that a consistent 
set of flags.
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
optics dependent. Tipically the user want to modify the as_built sequence by
adding special devices, markers, by rotating it or make it thin. All this
manipulation can be done there.
'''

# %%
ost.build_sequence(mad, beam=beam_to_configure)

# %%
'''
### Load a specific optics

This is the second example of user's defined function. Here the typically call the
stregth file.
'''

# %%
ost.apply_optics(mad, optics_file=optics_file)

# %%
'''
### Force disable beam-beam when needed
This as an example when the flag of the script superseed the mask_parameters.
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

!!! hint
	Having a link "module" in the working folder allow us to inspect the code in the 
	modules. This are in MAD-X code.
'''

# %%
mad.call("modules/submodule_01a_preparation.madx")
mad.call("modules/submodule_01b_beam.madx")

# %%
'''
### Check the machine of the repository

Here is important to note that the user's knob are not yet applied. 
Therefore the sequence(s) has(have) the knobs set from the optics repository.

To do this check the users can wrire their own function. Here you are an example.

!!! question
	TODO: Why do we put this function in the *ost*? 
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
### Set IP1/5 phase, apply and save crossing
'''

# %%
mad.call("modules/submodule_01c_phase.madx")

# %%
'''
### Set optics-specific knobs

Here the user needs to take care of the different conventions for the knob
definitions: some of the knobs definition are optics dependent (unfortunately). 
The user needs to make the proper knob conversion.  
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
### Check closed-orbit flatness
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
### Call the luminosity leveling module

In the following we show how to level the luminosity using a pythonic approach.
'''

def print_luminosity(mad, twiss_dfs, mask_parameters):
	for ip, number_of_ho_collisions in zip(['ip1', 'ip2', 'ip5', 'ip8'],
		[mask_parameters['par_nco_IP1'],
   	 mask_parameters['par_nco_IP2'],
   	 mask_parameters['par_nco_IP5'],
   	 mask_parameters['par_nco_IP8']]):
		myL=lumi.compute_luminosity(mad, twiss_dfs, ip, number_of_ho_collisions)
		print(f'L in {ip} is {myL} Hz/cm^2')

print_luminosity(mad, twiss_dfs, {k: mask_parameters[k] for k in ['par_nco_IP1',
																																	'par_nco_IP2',
																																	'par_nco_IP5',
																																	'par_nco_IP8']})

# %%
'''
#### Luminosity leveling by intensity in IP1/5. 

This is a trivial leveling, we use it just for the sake of the example.
'''

# %%
print('\n==== Intensity Luminosity Levelling ====')

from scipy.optimize import least_squares
L_target=2e+34
starting_guess=mad.sequence.lhcb1.beam.npart

def function_to_minimize(n_part):
	my_dict_IP1=lumi.get_luminosity_dict(mad, twiss_dfs,
 'ip1', mask_parameters['par_nco_IP1'])  
	my_dict_IP1['N1']=n_part
	my_dict_IP1['N2']=n_part
	my_dict_IP5=lumi.get_luminosity_dict(mad, twiss_dfs, 'ip5', mask_parameters['par_nco_IP5'])  
	my_dict_IP5['N1']=n_part
	my_dict_IP5['N2']=n_part
	return lumi.L(**my_dict_IP1)+lumi.L(**my_dict_IP5)-2*L_target

aux=least_squares(function_to_minimize, starting_guess)
print(aux)
print(f"\nLuminosity after levelling: {function_to_minimize(aux['x'][0])+L_target} Hz/cm^2")

mad.sequence.lhcb1.beam.npart=aux['x'][0]
mad.sequence.lhcb2.beam.npart=aux['x'][0]
mask_parameters['par_beam_npart']=aux['x'][0]

print_luminosity(mad, twiss_dfs, {k: mask_parameters[k] for k in ['par_nco_IP1',
																																	'par_nco_IP2',
																																	'par_nco_IP5',
																																	'par_nco_IP8']})

# %%
'''
#### Luminosity leveling by separation in IP8 

We separate the beams vertically in IP8.
'''

# %%
print('\n==== IP8 Luminosity Levelling ====')

L_target=mask_parameters['par_lumi_ip8']
# as starting guess it is good practice not to have a vanish derivative. 
# In fact, if the two beams are too much separated or if they are HO, 
# the algorithm could assume (wrongly) that the optimization is converging.
sigma_y_b1=np.sqrt(twiss_dfs['lhcb1'].loc['ip8:1'].bety*mad.sequence.lhcb1.beam.ey)
starting_guess=sigma_y_b1

def function_to_minimize(on_sep8v):
	my_dict_IP8=lumi.get_luminosity_dict(mad, twiss_dfs,
 'ip8', mask_parameters['par_nco_IP8'])  
	my_dict_IP8['y_1']=np.abs(on_sep8v)
	my_dict_IP8['y_2']=-np.abs(on_sep8v)
	return lumi.L(**my_dict_IP8)-L_target

aux=least_squares(function_to_minimize, starting_guess)
print(aux)
print(f"\nLuminosity after levelling: {function_to_minimize(aux['x'][0])+L_target} Hz/cm^2")

mad.globals['on_sep8v']=np.abs(aux['x'][0])*1e3

twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_after_ip8_leveling',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

mad.use(f'lhcb{beam_to_configure}')

mad.input('exec, crossing_save')
print('After IP8 leveling')
print_luminosity(mad, twiss_dfs, {k: mask_parameters[k] for k in ['par_nco_IP1',
																																	'par_nco_IP2',
																																	'par_nco_IP5',
																																	'par_nco_IP8']})

# %%
'''
#### Halo collision in IP2 

Here we consider a full horizontal separation of 5 sigmas in IP2.
'''

# %%
print('\n==== Halo collision in IP2 ====')

sigma_x_b1=np.sqrt(twiss_dfs['lhcb1'].loc['ip2:1'].betx*mad.sequence.lhcb1.beam.ex)
mad.globals['on_sep2h']=sigma_x_b1*5/2*1e3

twiss_dfs, other_data = ost.twiss_and_check(mad, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_after_ip2_halo_offset',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

mad.use(f'lhcb{beam_to_configure}')

mad.input('exec, crossing_save')
print('After IP2 halo offset')
print_luminosity(mad, twiss_dfs, {k: mask_parameters[k] for k in ['par_nco_IP1',
																																	'par_nco_IP2',
																																	'par_nco_IP5',
																																	'par_nco_IP8']})

# %%
'''
One could also use the old MAD-X approach.
'''

# %%
# if mode=='b4_without_bb':
#     print('Leveling not working in this mode!')
# else:
#     mad.call("modules/module_02_lumilevel.madx")

# %%
'''
Finally put at zero the dispersions bumps.
'''

# %%
mad.input('on_disp = 0')

# %%
'''
### Prepare the BB dataframes
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
	Here the dataframes can be edited, e.g., to set bbb intensity or to 
	follow the collision schedule of a specific bunch.
'''

# %%
'''
### Generate mad instance for B4

This is when the B2-B4 trasplant takes place.
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

From now on, we work exclusively on the sequence to track and we
refer to the *mad_track* handle for MAD-X.
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
### Twiss of the machine to track
'''

# %%
twiss_dfs, other_data = ost.twiss_and_check(mad_track, sequences_to_check,
        tol_beta=tol_beta, tol_sep=tol_sep,
        twiss_fname='twiss_track_intermediate',
        save_twiss_files=save_intermediate_twiss,
        check_betas_at_ips=check_betas_at_ips, check_separations_at_ips=False)

# %%
'''
### Installing the BB lenses

Here we install the BB lenses and depending on the flags we use the BB dataframes 
or the legacy approach (that is not working for B4).

The lenses are installed but their charge is off.
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

if enable_bb_legacy:
    assert(beam_to_configure == 1)
    assert(not(track_from_b4_mad_instance))
    assert(not(enable_bb_python))
    mad_track.call("modules/module_03_beambeam.madx")

# %%
'''
### Install crab cavities

Clearly this point is skipped for Run 3.
'''

# %%
# mad_track.call("tools/enable_crabcavities.madx")
# mad_track.call("modules/submodule_04_1a_install_crabs.madx")

# %%
'''
### Save references (orbit at IPs)

All orbits references are saved.
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
For the moment the errors of Run 3 are not implemented.
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

Nothing to do for Run 3.
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

# %%
'''
We hope that it is clear(er).
'''

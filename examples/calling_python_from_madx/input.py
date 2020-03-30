# %% Definition of the parameters that are not knobs of the beam sequence
parameter_dict={
    # =============================================================================
    # Beam parameters
    # =============================================================================
    ## LHC beam 1 (clockwise), LHC beam 2 (clockwise), LHC beam 2 (counterclockwise) 
    'par_mylhcbeam': 1, 
    ## beam normalized emittance [mm mrad]
    'par_beam_norm_emit': 2.5,
    ## [m]
    'par_beam_sigt': 0.075,
    ## [-]           
    'par_beam_sige': 1.1e-4,
    ## [-]                    
    'par_beam_npart': 1.4e11, 
    ## [GeV]            
    'par_beam_energy_tot': 7000,
    ## [A]          
    'par_oct_current': 350,
    ## [-]            
    'par_chromaticity': 15,
    ## [MV]          
    'par_vrf_total': 16.,
    ## Tunes with fractional part          
    'par_qx0': 62.31, 'par_qy0': 60.32,
    # =============================================================================
    # Beam-Beam configuration 
    # =============================================================================
    ## install the BB elements [0,1]
    'par_on_bb_switch': 1,
    ## if 1 lumi leveling in ip8 is applied and q/q' match is done with bb off [0,1]
    'par_on_collision': 1, 
    ## bunch separation [ns]               
    'par_b_t_dist': 25.,   
    ## default value for the number of additionnal parasitic encounters inside D1              
    'par_n_inside_D1': 5,                 
    ## number of slices for head-on in IR1 [between 0 and 201]
    'par_nho_IR1': 11, 'par_nho_IR2': 11, 'par_nho_IR5': 11, 'par_nho_IR8': 11, 
    ## flag to install the Crab Cavities [0, 1]
    'par_install_crabcavities': 0,
    # =============================================================================
    # Leveling in IP8   
    # =============================================================================
    # leveled luminosity in IP8 (considered if par_on_collision=1) [Hz/cm2]
    'par_lumi_ip8': 2e33,                 
    # These variables define the number of Head-On collisions in the 4 IPs
    'par_nco_IP1': 2592, 'par_nco_IP2': 2288, 'par_nco_IP5': 2592, 'par_nco_IP8': 2396,
    # =============================================================================
    # Errors and corrections 
    # =============================================================================
    # Select seed for errors
    'par_myseed': 0,
    # Set this flag to correct the errors of D2 in the NLC 
    # (warning: for now only correcting b3 of D2, still in development)
    'par_correct_for_D2': 0,
    # Set this flag to correct the errors of MCBXF in the NLC 
    # (warning: this might be less reproducable in reality, use with care)
    'par_correct_for_MCBX': 0,
    'par_off_all_errors': 0,
    'par_on_errors_LHC': 0,
    'par_on_errors_MBH': 0,
    'par_on_errors_Q5': 0,
    'par_on_errors_Q4': 0,
    'par_on_errors_D2': 0,
    'par_on_errors_D1': 0,
    'par_on_errors_IT': 0,
    'par_on_errors_MCBRD': 0,
    'par_on_errors_MCBXF': 0,
    # =============================================================================
    # Additional parameters
    # =============================================================================
    # parameter for having verbose output [0,1]
    'par_verbose': 1,
    # definition of the slicefactor used in the makethin
    'par_slicefactor': 4,
}

# %% Import few packages
from cpymad.madx import Madx
import numpy as np
import madxp
import os
from madxp import cpymadTool as mt

# %% Read the masked parameters
# In input.madx one can write the masked variable on the masked_variables.txt files
# Those are loaded on mask_dict
mask_dict={}
try:
    with open('masked_variables.txt', 'r') as parameter_file:
        for l in parameter_file:  
            parameter , value = l.strip().split('=')
            # forcing lower case for the parameter
            # forcing float for the values
            mask_dict[parameter.strip().lower()]=float(value.strip()) 
        print('File "masked_variables.txt" found.')
except FileNotFoundError:
    print('File "masked_variables.txt" not found. Custom data are assumed.')

# %% Overwrite if a variable, prensent in mask_dict, is present also in parameter_dict 
# the it is overwritten on parameter_dict. 
for i in mask_dict.keys():
    if i in parameter_dict.keys():
        print(f'Re-assign {i}')
        parameter_dict[i]=mask_dict[i]

# %%
command_log_file='log.madx'
stdout_file='stdout.madx'
with open(stdout_file, 'w') as myFile:
    mad = Madx(stdout=myFile,command_log=command_log_file)

# %%
working_folder='./'
os.chdir(working_folder)
my_df=madxp.madx2df('main.madx')
# %%

# Define integer tunes and tune split
parameter_dict['par_qx00']=int(parameter_dict['par_qx0'])
parameter_dict['par_qy00']=int(parameter_dict['par_qy0'])
parameter_dict['par_tsplit']=parameter_dict['par_qx00']-parameter_dict['par_qy00']

# to be checked on change of the parameter_dict
def check_parameter_dict(parameter_dict):
    assert parameter_dict['par_nco_IP5']==parameter_dict['par_nco_IP1']
    assert parameter_dict['par_qx00']-parameter_dict['par_qy00']==parameter_dict['par_tsplit']
    assert 'par_mylhcbeam' in parameter_dict
    assert 'par_beam_norm_emit' in parameter_dict

check_parameter_dict(parameter_dict)

for i in parameter_dict:
    mad.input(f'{i}={parameter_dict[i]};')
execution_df=madxp.df2run(mad,my_df['Prepare the environment':'Call check optics'])
# %% Some analysis 
# The execution DF has a lot of information
execution_df
# Let us focus on "Call check optics
# We have the dependentVariableDF
execution_df.loc['Call check optics'].dependent_variable_df
# We have the constant (less interesting)
execution_df.loc['Call check optics'].constant_df
# We the independent variables
execution_df.loc['Call check optics'].independent_variable_df
# We have the table list at that moment
execution_df.loc['Call check optics'].tables_list
# that we can access with (IF not overwritten)
mt.summ_df(mad.table.summ)
mt.twiss_df(mad.table.twiss_bare_optics_b1)
# The sequences
execution_df.loc['Call check optics'].sequences_df
# We can get DF of a sequence (IF not overwritten, no drift present
my_sequence_df=mt.sequence_df(mad, 'lhcb1')
my_sequence_df
# and finally one can have the knobs
mt.knobs_df(my_sequence_df)
# check all the element of a sequence depending on a given knob
mt.knob_df('on_x1',my_sequence_df)
# on can also use it for the dependentVariableDF
mt.knob_df('on_x1',execution_df.loc['Call check optics'].dependent_variable_df)
# and one can analyze a 
mt.show_element('ms.11r3.b1',my_sequence_df)
# The beams attached for the sequences
execution_df.loc['Call check optics'].beams_df

# %% Saving and reading in pandas
import pandas as pd
mt.table_df(mad.table['twiss_bare_optics_b1']).to_parquet(f'twiss_bare_optics_b1_{parameter_dict["par_slicefactor"]}.parquet')
display(pd.read_parquet(f'twiss_bare_optics_b1_{parameter_dict["par_slicefactor"]}.parquet',columns=['betx','bety']))

# %% Simple checks
ip_list=['ip1:1','ip2:1','ip5:1','ip8:1']
columns_list=['s','betx','alfx','bety','alfy']
twiss_df_b1=mt.twiss_df(mad.table.twiss_bare_optics_b1)
display(twiss_df_b1.loc[ip_list][columns_list])
twiss_df_b2=mt.twiss_df(mad.table.twiss_bare_optics_b2)
display(twiss_df_b2.loc[ip_list][columns_list])

# %% Make some plots
from matplotlib import pylab as plt
for my_table, my_title in zip(['twiss_bare_optics_b1',
                            'twiss_bare_optics_b2'],
                            ['B1, from optics repository', 
                            'B2, from optics repository']):
    aux = mt.table_df(mad.table[my_table])
    plt.figure(figsize=(15,10))
    plt.plot(aux.s,aux.x,'b', lw=1, label='x')
    plt.plot(aux.s,aux.y,'r',lw=1, label='y')
    plt.title(my_title)
    plt.xticks([aux.loc['ip3:1'].s,
        aux.loc['ip4:1'].s,
        aux.loc['ip5:1'].s,
        aux.loc['ip6:1'].s,
        aux.loc['ip7:1'].s,
        aux.loc['ip8:1'].s,
        aux.loc['ip1:1'].s,
        aux.loc['ip2:1'].s],['IP'+str(x) for x in [3,4,5,6,7,8,1,2]])
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylabel('[m]')

# %%
def from_beta_to_xing_angle_urad(beta_m):
    return  0.5*(132.47 + 58.3959 * np.sqrt(beta_m) + 30.0211 * beta_m)/np.sqrt(beta_m)

knob_dict={
    'on_sep1': 0,  
    'on_sep5': 0,         
    'on_sep2h': 0,
    'on_sep2v': 0,
    'on_x2h': 0,
    'on_x2v': 200,
    'on_sep8h': 0,
    'on_sep8v': 0,
    'on_x8h': 0,
    'on_x8v': 135,
    'par_on_qpp': 0,
    'on_disp': 1,
    'on_alice': 7000/parameter_dict['par_beam_energy_tot'],
    'on_lhcb': 7000/parameter_dict['par_beam_energy_tot'],
    'on_sol_atlas': 7000/parameter_dict['par_beam_energy_tot'],
    'on_sol_cms': 7000/parameter_dict['par_beam_energy_tot'],
    'on_sol_alice': 7000/parameter_dict['par_beam_energy_tot'],
}
betx_ip1 = mad.globals['betx_ip1']
knob_dict['on_x1'] = from_beta_to_xing_angle_urad(betx_ip1)
knob_dict['on_x5'] = from_beta_to_xing_angle_urad(betx_ip1)

for i in mask_dict.keys():
    if i in knob_dict.keys():
        print(f'Set {i}')
        knob_dict[i] = mask_dict[i]

for i in knob_dict:
    mad.input(f'{i} = {knob_dict[i]};')

#    !Avoid crabbing more than the crossing angle
#    if ( abs(on_crab1)>abs(par_xing_ang_ip15) && on_crab1 <> 0) {on_crab1 = abs(on_crab1)/on_crab1 * abs(par_xing_ang_ip15);}
#    if ( abs(on_crab5)>abs(par_xing_ang_ip15) && on_crab5 <> 0) {on_crab5 = abs(on_crab5)/on_crab5 * abs(par_xing_ang_ip15);}

# %%
import pandas as pd
aux_df = madxp.df2run(mad,my_df['Call define and save crossing':'Call define and save crossing'])
execution_df = pd.concat([execution_df, aux_df])

# %% Usual plots on orbits
for my_table, my_title in zip(['twiss_crossing_disable_b1',
                            'twiss_crossing_disable_b2',
                            'twiss_crossing_enable_b1',
                            'twiss_crossing_enable_b2'],
                            ['B1, crossing disabled', 
                            'B2, crossing disabled',
                            'B1, crossing enabled',
                            'B2, crossing enabled']):
    aux = mt.table_df(mad.table[my_table])
    plt.figure(figsize=(15,10))
    plt.plot(aux.s,aux.x,'b', lw=1, label='x')
    plt.plot(aux.s,aux.y,'r',lw=1, label='y')
    plt.title(my_title)
    plt.xticks([aux.loc['ip3:1'].s,
        aux.loc['ip4:1'].s,
        aux.loc['ip5:1'].s,
        aux.loc['ip6:1'].s,
        aux.loc['ip7:1'].s,
        aux.loc['ip8:1'].s,
        aux.loc['ip1:1'].s,
        aux.loc['ip2:1'].s],['IP'+str(x) for x in [3,4,5,6,7,8,1,2]])
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylabel('[m]')

# %% Usual plots on optics
for my_table, my_title in zip(['twiss_crossing_enable_b1',
                            'twiss_crossing_enable_b2'],
                            ['B1, crossing enabled',
                            'B2, crossing enabled']):
    aux = mt.table_df(mad.table[my_table])
    aux = aux.loc['ip4:1':'ip6:1'] # smart filtering
    plt.figure(figsize=(10,7))
    plt.plot(aux.s,aux.betx,'b', lw=1, label='betx')
    plt.plot(aux.s,aux.bety,'r',lw=1, label='bety')
    plt.title(my_title)
    plt.xticks([aux.loc['ip4:1'].s,
            aux.loc['ip5:1'].s,
            aux.loc['ip6:1'].s],
            ['IP'+str(x) for x in [4,5,6]])
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylabel('[m]')
# %% Plotting the bending and gradients
for beam in ['1','2']:
    aux = mt.table_df(mad.table['twiss_crossing_enable_b'+beam])
    aux=aux.loc['e.ds.l1.b'+beam+':1':'s.ds.r1.b'+beam+':1']
    plt.figure(figsize=(20,3))
    plt.fill_between(aux['s'],aux['k0l']*1e3,step="pre",color='b')
    ax1=plt.gca()
    ax1.set_ylabel('k0l [mrad]', color='b') 
    ax1.tick_params(axis='y', labelcolor='b')
    ax2 = ax1.twinx() 
    ax2.fill_between(aux['s'],aux['k1l']*1e3,step="pre",color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.set_ylabel('k1l [km$^{-1}$]', color='r') 
    plt.title('BEAM '+beam)
    ax1.set_xlabel('s from IP3 [m]')

# %% Plotting the bending and gradients
mad.input('use,sequence=lhcb1; survey;')
# %%

plt.plot(mt.table_df(mad.table.survey_b1))

# %%
aux_df = madxp.df2run(mad,my_df['Luminosity levelling':'Luminosity levelling'])
execution_df = pd.concat([execution_df, aux_df])




# %% From Sofia
def cmp_luminosity_0(intensity, bunches, sigma_x, sigma_P, frev=11245.5):
    return (intensity**2*bunches*frev)/(4.0*np.pi*sigma_x*sigma_P)

def crossing_factor(phi, sigmaz, sigma_x):
    return 1./np.sqrt(1. + (sigmaz*phi/2.)**2/(sigma_x**2)  )

def offset_factor(halo):
    a  = -(halo/ 2.)**2
    return np.exp(a)

def cmp_luminosity(intensity, bunches, sigma_x, sigma_P, phi, sigmaz, halo,frev=11245.5):
    crossing = crossing_factor(phi, sigmaz, sigma_x)
    offset   = offset_factor(halo)
    #print 'Initial lumi is : ', cmp_luminosity_0(intensity, bunches, sigma_x, sigma_P, frev)*1e-4
    #print 'Crossing factor is: ',crossing
    #print 'Offset factor is: ',offset
    return cmp_luminosity_0(intensity, bunches, sigma_x, sigma_P, frev)*crossing*offset

intensity_b1 = mad.sequence['lhcb1'].beam['npart']
bunch_length_b1 = mad.sequence['lhcb1'].beam['sigt']
intensity_b2 = mad.sequence['lhcb2'].beam['npart']
bunch_length_b2 = mad.sequence['lhcb2'].beam['sigt']

assert  intensity_b1 == intensity_b2
assert  bunch_length_b2 == bunch_length_b2

def compute_luminosity(nIP, mad, twiss_DF_b1, twiss_DF_b2):
    bunches = mad.globals['par_nco_IP8']

    myDict={}
    myDict['b1']={}
    myDict['b1']['DF']=twiss_DF_b1.loc[[f'ip{nIP}:1']]
    myDict['b2']={}
    myDict['b2']['DF']=twiss_DF_b2.loc[[f'ip{nIP}:1']]


    myDict['b1']['ex']  = mad.sequence['lhcb1'].beam['ex']
    myDict['b1']['ey']  = mad.sequence['lhcb1'].beam['ey']
    myDict['b1']['betx']= myDict['b1']['DF'].loc[f'ip{nIP}:1']
    myDict['b1']['bety']= myDict['b1']['DF'].loc[f'ip{nIP}:1']

    myDict['b1']['sigmax'] = np.sqrt(myDict['b1']['ex']*myDict['b1']['DF'].betx)
    myDict['b1']['sigmay'] = np.sqrt(myDict['b1']['ey']*myDict['b1']['DF'].betx)
    cmp_luminosity_0(intensity_b1, bunches, sigma_x, sigma_P, frev=11245.5)
    return myDict


compute_luminosity(8, mad, mt.table_df(mad.table.twiss_crossing_enable_b1), mt.table_df(mad.table.twiss_crossing_enable_b2))
# %%
mad.input('quit;')
print('Exit python.')
# %%

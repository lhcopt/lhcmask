# %% Import packages
from cpymad.madx import Madx

# %% Check the executable
import sys
print(sys.executable)

# %% Changing folder
import os
# please change it conveniently

# os.chdir('/afs/cern.ch/work/s/sterbini/beambeam_macros/examples')
print(os.getcwd())
# %% Unmask the mask
fname_mask = 'hl14_collisions.mask'

def unmask_file(fname_mask,fname_unmasked,
beam = 1,
i_octupoles = 100,
emittance_in_m = 2.3e-6,
n_particles = 2.25e11,
chromaticity = 15,
xing_angle_urad = 245,
seedran = 1,
):

    with open(fname_mask) as fid:
        mask_content = fid.read()

    mask_content = mask_content.replace(r'%BEAM%', str(1))
    mask_content = mask_content.replace(r'%OCT%', f'{i_octupoles:e}')
    mask_content = mask_content.replace(r'%EMIT_BEAM', f'{emittance_in_m:e}')
    mask_content = mask_content.replace(r'%NPART', f'{n_particles:e}')
    mask_content = mask_content.replace(r'%CHROM%', f'{chromaticity:e}')
    mask_content = mask_content.replace(r'%XING', f'{xing_angle_urad:e}')
    mask_content = mask_content.replace(r'%SEEDRAN', f'{seedran:d}')
    with open(fname_unmasked, 'w') as fid:
        fid.write(mask_content)
# %% Dump the unmasked mask on file
#unmask_fname = fname_mask.split('.mask')[0]+'_unmask.mask'
path='./collision/'
unmask_file(fname_mask=path+'new/hl14_collisions.mask',fname_unmasked=path+'new/hl14_collisions.madx')
unmask_file(fname_mask=path+'old/old.mask',fname_unmasked=path+'old/old.madx')

# %%
def f_proc1():
    mad = Madx()
    mad.chdir(path+'old')
    mad.call('old.madx')

# %%
def f_proc2():
    mad = Madx()
    mad.chdir(path+'new')
    mad.call('hl14_collisions.madx')

from multiprocessing import Process

p1 = Process(target=f_proc1)
p2 = Process(target=f_proc2)

p1.start()
p2.start()

p1.join()
p2.join()

print('Mah!')

# %%

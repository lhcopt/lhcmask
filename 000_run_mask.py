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
beam = 1
i_octupoles = 100.
emittance_in_um = 2.3
n_particles = 2.25e11
chromaticity = 15
xing_angle_urad = 245.
seedran = 1

fname_mask = 'hl14_collisions.mask'

with open(fname_mask) as fid:
    mask_content = fid.read()

mask_content = mask_content.replace(r'%BEAM%', str(1))
mask_content = mask_content.replace(r'%OCT%', f'{i_octupoles:e}')
mask_content = mask_content.replace(r'%EMIT_BEAM', f'{emittance_in_um:e}')
mask_content = mask_content.replace(r'%NPART', f'{n_particles:e}')
mask_content = mask_content.replace(r'%CHROM%', f'{chromaticity:e}')
mask_content = mask_content.replace(r'%XING', f'{xing_angle_urad:e}')
mask_content = mask_content.replace(r'%SEEDRAN', f'{seedran:d}')
# %% Dump the unmasked mask on file
unmask_fname = fname_mask.split('.mask')[0]+'_unmask.mask'
with open(unmask_fname, 'w') as fid:
    fid.write(mask_content)

#os.system('madx '+unmask_fname)
mad = Madx()
mad.call(unmask_fname)
# %%

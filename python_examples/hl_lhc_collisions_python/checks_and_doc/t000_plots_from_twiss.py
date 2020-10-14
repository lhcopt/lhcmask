import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

twiss_df = pd.read_parquet('../twiss_with_crossing_seq_lhcb1.parquet')


plt.close('all')

# %%
fig = plt.figure(1, figsize=(6.4*1.6, 4.8*1.5))
axbetx = fig.add_subplot(2,2,1)
axbetx.plot(twiss_df['s'], twiss_df['betx'])
axbetx.set_xlabel('s')
axbetx.set_ylabel(r'$\beta_x$ [m]')
axbetx.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')

# %%
axbety = fig.add_subplot(2,2,2, sharex=axbetx, sharey=axbetx)
axbety.plot(twiss_df['s'], twiss_df['bety'])
axbety.set_xlabel('s')
axbety.set_ylabel(r'$\beta_y$ [m]')
axbety.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
axbety.legend()
# %%

axcox = fig.add_subplot(2,2,3, sharex=axbetx)
axcox.plot(twiss_df['s'], twiss_df['x'])
axcox.set_xlabel('s')
axcox.set_ylabel('x [m]')
#plt.xlim(5500,7500)
# %%
axcoy = fig.add_subplot(2,2,4, sharex=axbetx, sharey=axcox)
axcoy.plot(twiss_df['s'], twiss_df['y'])
axcoy.set_xlabel('s')
axcoy.set_ylabel('y [m]')

plt.show()

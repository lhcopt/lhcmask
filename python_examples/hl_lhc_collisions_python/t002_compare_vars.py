import pickle

fname1 = './vars_end.pkl'
fname2 = '../hl_lhc_collision_nobb_b4/vars_end.pkl'

with open(fname1, 'rb') as fid:
    d1 = pickle.load(fid)
with open(fname2, 'rb') as fid:
    d2 = pickle.load(fid)

group = 'dependent_variables_expr'
# group = 'constants'
# group = 'independent_variables'
for kk in d1[group].keys():
    if kk in d2[group].keys():
        if d1[group][kk] != d2[group][kk]:
            print(kk, d1[group][kk], d2[group][kk])
print('FOR MAD-X:')
for kk in d1[group].keys():
    if kk in d2[group].keys():
        if d1[group][kk] != d2[group][kk]:
            print(f'{kk} = {d1[group][kk]:.50e};')

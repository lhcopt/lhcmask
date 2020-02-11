
def unmask(mask_filename, parameters,
        output_filename=None, strict=False):

    with open(mask_filename, 'r') as fid:
        content = fid.read()

    for kk in parameters.keys():
        content = content.replace(kk, str(parameters[kk]))

    if output_filename is None:
        outfname = None
    elif output_filename == 'auto':
        outfname = mask_filename+'.unmasked'
    else:
        outfname = output_filename

    if outfname is not None:
        with open(outfname, 'w') as fid:
            fid.write(content)

    if strict:
        if '%%' in content:
            raise ValueError(
              ('There is still a %% after unmasking!\n'
               '--> Incompatible with strict=True'))

    return content

def parse_parameter_file(fname):
    with open(fname, 'r') as fid:
        lines = fid.readlines()

    ddd = {}
    for ii, ll in enumerate(lines):
        if ':' not in ll:
            print('Warning! line %d is skipped! Its content is:'%(ii+1))
            print(ll)
            continue

        nn = ll.split(':')[0].replace(' ', '').replace('\n', '')
        vv = ll.split(':')[1].replace(' ', '').replace('\n', '')
        ddd[nn] = vv

    return ddd

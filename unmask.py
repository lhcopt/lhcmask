
def unmask(mask_filename, parameters,
        output_filename=None):

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

    return content


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Unmasks a madx maskfile.')
    parser.add_argument('mask_filename', help='Name of the maskfile to be unmasked')
    parser.add_argument('parameters', nargs='+', help=('Parameters for unmasking (can be a filename,'
                    'or given inline as: %%PARAM1%%:value1, %%PARAM2:value2, ...'))
    parser.add_argument('--output_filename', help='Name of the output file', default='auto')
    args = parser.parse_args()

    if ':' in args.parameters[0]:
        par_dict = {}
        for pp in args.parameters:
            kk = pp.split(':')[0].replace(' ', '')
            vv = pp.split(':')[1].replace(' ', '')
            par_dict[kk] = vv
    else:
        pass

    unmask(args.mask_filename, par_dict, output_filename=args.output_filename)

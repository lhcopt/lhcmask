
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

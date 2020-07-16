import os 

if __name__=='__main__':
	print('Starting py to ipynb conversion...', end='')
	import argparse
	parser = argparse.ArgumentParser(description='From a py file make a ipynb and a md file.')
	parser.add_argument('mask_filename', help='Name of the maskfile to be unmasked')
	args = parser.parse_args()
	filename=args.mask_filename.split('.')[0]
	os.system(f'ipynb-py-convert {filename}.py {filename}.ipynb')
	print('done.')
	print('Starting ipynb to md conversion.', flush=True)
	os.system(f'jupyter nbconvert --to markdown {filename}.ipynb')
	print('Done.')

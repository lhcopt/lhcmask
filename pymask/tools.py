import ruamel.yaml
import datetime

# load the configuration
def read_yaml(myfile, verbose=True):
	"""
	read_yaml is to read a yaml file and convert it into python.
	Example
	--------
	>>> read_yaml('mytest')
	FileNotFoundError, empty dict
	{}
	"""
	try: 
		with open(myfile, 'r') as file:
			yaml = ruamel.yaml.YAML()    
			my_dict = yaml.load(file)
		return my_dict
	except FileNotFoundError:
		if verbose: print('FileNotFoundError, empty dict')
		my_dict = {}
		return my_dict
	except ruamel.yaml.constructor.DuplicateKeyError:
		my_dict = {}
		return my_dict
	except Exception as e: 
		print(e.__class__)
		return None


def write_yaml(my_dict, myfile):
	"""
	write_yaml is to convert a dictionary into a yaml file.
	
	Examples
	--------
	>>> write_yaml({'green': 'hello'}, ('mytest2.yml'))
	"""	
	with open(myfile, 'w') as file:  
		yaml = ruamel.yaml.YAML()     
		yaml.dump(my_dict, file)


def append_yaml(my_dict, myfile):
	"""
	append_yaml is to append dictionaries to a yaml file.
	Examples
	--------
	>>> append_yaml({'blue': 'bonjour'}, ('mytest2.yml'))
	"""	
	with open(myfile, 'a') as file:  
		yaml = ruamel.yaml.YAML()    
		yaml.dump(my_dict, file)



def get_last_stage(myfile, verbose=True):
	"""
	get_last_stage
	Examples
	--------
	>>> tag_it('myfile', 'hello')
	"""	
	my_dict=read_yaml(myfile, verbose)
	try:
		return list(my_dict.keys())[-1]+1
	except IndexError:
		if verbose: print('IndexError, I consider 0 as first item')
		return 0
	except Exception as e: 
		print(e.__class__)
		return 0


def tag_it(myfile, mycomment):
	"""
	tag_it is to create timestamps and add them to a yaml file.
	Examples
	--------
	>>> tag_it('myfile', 'hello')
	"""	
	stage = get_last_stage(myfile)
	with open(myfile, 'a') as file:
		yaml = ruamel.yaml.YAML() 
		my_dict = {stage: {}}
		my_dict[stage]['human_time'] = datetime.datetime.now()
		my_dict[stage]['unix_time'] = datetime.datetime.now().timestamp()        #in seconds
		my_dict[stage]['comment'] = mycomment
		yaml.dump(my_dict, file)
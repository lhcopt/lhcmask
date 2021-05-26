import ruamel.yaml
import datetime

# load the configuration
def read_yaml(myfile, verbose=True):
	"""
	read_yaml is to read a yaml file and convert it into python.
	Example
	--------
	>>> read_yaml('mytest.yml')

	This will return the contents inside of a yaml file.
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

	This will add {'green': 'hello'} as a dictionary to the ('mytest2.yml') file.
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

	This will append {'blue': 'bonjour'} to an existing yaml file: ('mytest2.yml')
	"""	
	with open(myfile, 'a') as file:  
		yaml = ruamel.yaml.YAML()    
		yaml.dump(my_dict, file)



def get_last_stage(myfile, verbose=True):
	"""
	get_last_stage is to read a yaml file and to return the number of the last dictionary

	Examples
	--------
	>>> get_last_stage('myfile', verbose=True)

	If the get_last_stage has two dictionaries: labeled '0' and '1' it will return the '1', the last stage.
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

	This will create a a human readable and a unix time stamp and append that to 'myfile', includin the comment 'hello'
	"""	
	stage = get_last_stage(myfile)
	with open(myfile, 'a') as file:
		yaml = ruamel.yaml.YAML() 
		my_dict = {stage: {}}
		my_dict[stage]['human_time'] = datetime.datetime.now()
		my_dict[stage]['unix_time'] = datetime.datetime.now().timestamp()        #in seconds
		my_dict[stage]['comment'] = mycomment
		yaml.dump(my_dict, file)


def tag_first(myfile, mycomment):
	"""
	tag_first is to overwrite an already existing yaml file or create a new yaml file and add the first timestamps.

	Examples
	--------
	>>> tag_it('myfile', 'hello')

	If 'my_file' already exits, it will overwrite that file and create a first timestamp.
	"""	
	stage = 0
	with open(myfile, 'w') as file:
		yaml = ruamel.yaml.YAML() 
		my_dict = {stage: {}}
		my_dict[stage]['human_time'] = datetime.datetime.now()
		my_dict[stage]['unix_time'] = datetime.datetime.now().timestamp()        #in seconds
		my_dict[stage]['comment'] = mycomment
		yaml.dump(my_dict, file)


"""Installation with setuptools or pip."""
from setuptools import setup, find_packages
import os
import ast


def get_version_from_init():
    """Obtain library version from main init."""
    init_file = os.path.join(
        os.path.dirname(__file__), 'pymask', '__init__.py'
    )
    with open(init_file) as fd:
        for line in fd:
            if line.startswith('__version__'):
                return ast.literal_eval(line.split('=', 1)[1].strip())


setup(
    name='pymask',
    version=get_version_from_init(),
    description='Mask modules for LHC and HL-LHC',
    author='Guido Sterbini',
    author_email='guido.sterbini@cern.ch',
    url='https://github.com/lhcopt/lhcmask',
    license=None,
    packages=find_packages(exclude=('tests', 'docs', 'examples')),
    install_requires=[
        'numpy',
        'cpymad',
        'pandas',
    ],
)

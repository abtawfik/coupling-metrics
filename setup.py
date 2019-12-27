from setuptools import setup  #, find_packages


entry = '''
[console_scripts]
comet=comet.cli:comet
'''

setup( name = 'CoMeT',
       version='1.0',
       author_email='abtawfik@umich.edu',
       packages=['comet'],
       include_package_data=True,
       install_requires=['click', 'xarray', 'dask[complete]', 'toolz'],
       entry_points=entry
     )

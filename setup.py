from setuptools import setup


entry = '''
[console_scripts]
comet=comet.cli:comet
'''

setup( name = 'CoMeT',
       version='1.0',
       author_email='abtawfik@umich.edu',
       packages=['comet', 'comet.metrics'],
       include_package_data=True,
       install_requires=['click',
                         'xarray',
                         'dask[complete]',
                         'toolz',
                         'numpy',
                         'netcdf4',
                         'cfgrib'],
       entry_points=entry
     )

from setuptools import setup
from numpy.distutils.core import setup, Extension



lib = [Extension(name='comet.metrics.fortran.hcf'   , sources=['comet/metrics/fortran/hcfcalc.f90']),
       Extension(name='comet.metrics.fortran.ctp'   , sources=['comet/metrics/fortran/ctp_hilow.f90']),
       Extension(name='comet.metrics.fortran.rhtend', sources=['comet/metrics/fortran/rh_tendency.f90'])]


entry = '''
[console_scripts]
comet=comet.cli:comet
'''

setup( name = 'CoMeT',
       version='1.0',
       author='Ahmed Tawfik',
       author_email='abtawfik@umich.edu',
       url='https://github.com/abtawfik/coupling-metrics',
       packages=['comet', 'comet.metrics'],
       ext_modules=lib,
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

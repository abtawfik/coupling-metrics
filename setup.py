#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages, setup

#########################################
# Function for reading long description #
#########################################
def read(*names, **kwargs):
    with io.open( join(dirname(__file__), *names),
                  encoding=kwargs.get('encoding', 'utf8')
                 ) as fh:
        return fh.read()


####################
# Long Description #
####################
long_desc = '%s\n%s' % (re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')), re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst')))

################
# Project URLs #
################
project_urls = {'Documentation': 'https://coupling-metrics.readthedocs.io/',
                'Changelog': 'https://coupling-metrics.readthedocs.io/en/latest/changelog.html',
                'Issue Tracker': 'https://github.com/abtawfik/coupling-metrics/issues'}

#####################
# required packages #
#####################
install_reqs = ['click',
                'xarray',
                'dask[complete]',
                'toolz',
                'numpy',
                'netcdf4',
                'cfgrib']



#################
# Setup package #
#################
setup(name='coupling-metrics',
      version='1.0.0',
      license='MIT License',
      description='Calculate various state-of-the-art land-atmosphere coupling metrics',
      long_description=long_desc,
      author='Ahmed Tawfik',
      author_email='abtawfik@umich.edu',
      url='https://github.com/abtawfik/coupling-metrics',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
      include_package_data=True,
      zip_safe=False,
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Developers',
                   'Operating System :: Unix',
                   'Operating System :: POSIX',
                   'Operating System :: Microsoft :: Windows',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: Implementation :: CPython',
                   'Programming Language :: Python :: Implementation :: PyPy'],
      project_urls=project_urls,
      python_requires='>=3.*',
      install_requires=install_reqs,
      entry_points={'console_scripts': ['comet = comet.cli:comet']})



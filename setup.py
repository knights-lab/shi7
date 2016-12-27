from setuptools import setup, find_packages
from glob import glob
import os
import platform

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Kaiwei Ang", "Gabe Al-Ghalith", "Benjamin Hillmann"]
__email__ = "hillmannben@gmail.com"
__license__ = "MIT"
__maintainer__ = "Kaiwei Ang"
__version__ = "0.0.1-dev"

long_description = 'All of your shi7 is one place. From quality scores to mappable reads.'

if platform.system() == 'Darwin':
    print('Darwin found!')
else:
    print('Darwin not found!')

setup(
    name='shi7en',
    version=__version__,
    packages=find_packages(),
    package_data = {'shi7en': ['adapters/*.fa']},
    url='',
    license=__license__,
    author=__author__,
    author_email=__email__,
    description='',
    long_description=long_description,
    keywords='',
    install_requires=[],
    entry_points={
        'console_scripts': [
            'shi7en = shi7en.shi7en:main',
        ]
    },
)

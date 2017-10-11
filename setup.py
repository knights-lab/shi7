from setuptools import setup, find_packages
from glob import glob
import os
import platform
import versioneer

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Kaiwei Ang", "Gabe Al-Ghalith", "Benjamin Hillmann"]
__email__ = "hillmannben@gmail.com"
__license__ = "AGPL"
__maintainer__ = "Benjamin Hillmann"

long_description = 'All of your shi7 in one place. From quality scores to mappable reads.'

setup(
    name='shi7',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'shi7': ['adapters/*.fa']},
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
            'shi7 = shi7.shi7:main',
            'shi7_learning = shi7.shi7_learning:main'
        ]
    },
)

# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from setuptools import setup, find_packages
import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join(path, filename)[22:])
    return paths

extra_files = package_files(os.path.join('.', 'turbulucid', 'datasets'))

setup(name='turbulucid',
      version='0.2',
      description='A package for post-processing of 2D flow fields.',
      url='https://github.com/timofeymukha/turbulucid',

      author='Timofey Mukha',
      author_email='timofey.mukha@it.uu.se',
      packages=find_packages(),
      scripts=[],
      entry_points={
          'console_scripts': [
              'averageAlongAxis=turbulucid.bin.averageAlongAxis:main'
          ]
      },
      include_package_data=True,
      package_data={
          'turbulucid.datasets':extra_files
      },
      install_requires=[
                        'numpy',
                        'scipy',
                        'matplotlib',
                        'pytest'
                       ],
      license="GNU GPL 3",
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
      ],
      zip_safe=False)


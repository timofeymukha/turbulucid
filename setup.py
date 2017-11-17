# This file is part of turbulucid
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

from setuptools import setup

setup(name='turbulucid',
      version='0.1',
      description='A package for post-processing of 2D flow fields.',
      url='https://github.com/timofeymukha/turbulucid',

      author='Timofey Mukha',
      author_email='timofey.mukha@it.uu.se',
      packages=['turbulucid'],
      scripts=[],
      entry_points={
          'console_scripts': [
              'averageAlongAxis=turbulucid.bin.averageAlongAxis:main'
          ]
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


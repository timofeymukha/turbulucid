from setuptools import setup

setup(name='turbulucid',
      version='0.1',
      description='A package for simple post-processing of 2d turbulent fluid fields.',
      url='http://://bitbucket.org/lesituu/turbulucid',
      author='Timofey Mukha',
      author_email='timofey.mukha@it.uu.se',
      packages=['turbulucid'],
      scripts=[],
      install_requires=[
                        'numpy',
                        'scipy',
                        'matplotlib',
                        'pandas',
                        'PyFoam',
                        'h5py'
                       ],
      zip_safe=False)


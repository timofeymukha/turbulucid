from setuptools import setup

setup(name='turbulucid',
      version='0.1',
      description='A package for simple post-processing 2d turbulent fluid fields.',
      url='http://://bitbucket.org/lesituu/turbulucid',
      author='Timofey Mukha',
      author_email='TimofeyMukha@it.uu.se',
      packages=['turbulucid'],
      scripts=[],
      install_requires=[
                    'numpy',
                    'matplotlib'
                       ],
      zip_safe=False)


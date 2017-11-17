User Guide
==========

Prerequisites
-------------
Turbulucid is mainly tested using Python 3, but is known to work with Python 2
as well.
It should work on any platform, where the packages described below also work.
It has been tested extensively on both Windows and Linux.

For turbulucid to work, several other packages have to installed.
Three packages are :code:`numpy`, :code:`scipy` and :code:`matplotlib`.
These are easy to obtain and are part of many python distributions, in
particular, Anaconda.

Turbulucid also depends on python bindings for VTK, i.e the :code:`vtk`
 package.
The version of VTK should be at least 7.0.0.
Hopefully later versions work as well, but this has not been tested.
With Anaconda VTK can be obtained by running
:code:`conda install -c conda-forge vtk` in the terminal.

Installing
----------
Installing the package is easy.
Simply clone the git repository or download it as an archive and then
unpack.
Then navigate to the root catalog of the code in a terminal and execute
:code:`python setup.py install`.














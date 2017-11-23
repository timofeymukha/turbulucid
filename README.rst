turbulucid
==========

.. image:: https://travis-ci.org/timofeymukha/turbulucid.svg?branch=master
    :target: https://travis-ci.org/timofeymukha/turbulucid

What?
-----

Turbulucid is a package for post-processing two-dimensional cell-centered VTK
polyData.
The main use case is envisioned to be analysis of cut-plane data coming from
finite-volume based solvers for Computational Fluid Dynamics.
The package contains functions for both plotting the data and data
introspection.

Why?
----

VTK has become a popular format for storing unstructured datasets
(both 3d and 2d).
The API of VTK and software such as Paraview already provide the means for
working with VTK data, so why the need for a new package?
The answer is that while the above-mentioned tools are excellent for general
inspection of large 3d datasets, they do not provide the means for producing
*publication-quality* plots of 2d data.
An example of such a plot, produced with turbulucid, is shown below.

.. _fig-cover:

.. figure:: docs/figures/cover.png
   :align: center

Also, while VTK provides an abundance of filters for extracting specific
parts of a dataset, the object-oriented API is hard to learn and quite verbose.
Turbulucid provides easier access to the data and a set of functions for
performing simple extractions.

How?
----

Under the hood turbulucid uses VTK to handle the data, but exposes everything
to the user in terms of numpy arrays.
The plotting functions use matplotlib, and return associated matplotlib
objects.
This allows the user to harness the full customization power of matplotlib
to make the plots look exactly as desired.

Installing
----------
Turbulucid is mainly tested using both Python 3 and 2, but Python 3
is recommended.
It should work on any platform, where the packages described below also work.
It has been used extensively on both Windows and Linux.

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

Installing the package is easy.
Simply clone the git repository or download it as an archive and then
unpack.
Then navigate to the root catalog of the code in a terminal and execute
:code:`python setup.py install`.

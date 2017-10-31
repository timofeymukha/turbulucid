turbulucid
==========

Turbulucid is a package for post-processing two-dimensional CFD data.
The should be saved in VTK format, and be cell-centred.
This latter is typical of data produced by finite-volume based solvers.

Turbulucid provides functionality for both performing computations on the
data as well as producing high-quality plots.
This is achieved by providing simple to use data-extraction functions, that
use VTK as a back-bone, and plotting functions that use a combination of
VTK and matplotlib.

.. toctree::
   :caption: Table of Contents

   code_reference/index.rst

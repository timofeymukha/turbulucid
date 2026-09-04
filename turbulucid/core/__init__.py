# This file is part of turbulucid
# (c) 2018 Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the README for more information

from . import case, data_extraction, plotting, quantities, readers
from .case import Case
from .data_extraction import (
    dist,
    edge_lengths,
    isoline,
    normals,
    profile_along_line,
    sample_by_plane,
    sort_indices,
    tangents,
)
from .plotting import (
    add_colorbar,
    plot_boundaries,
    plot_contour,
    plot_field,
    plot_streamlines,
    plot_vectors,
)
from .quantities import delta_99, delta_star, momentum_thickness
from .readers import (
    LegacyReader,
    NativeReader,
    Reader,
    VTUReader,
    XMLReader,
    mark_boundary_cells,
)

# Kept as a literal so that static tooling can see it. The submodule
# order matches the imports above, and test_public_api guards it against
# drifting away from the submodules' own __all__.
__all__ = [
    "case", "plotting", "data_extraction", "quantities", "readers",
    # case
    "Case",
    # plotting
    "plot_boundaries", "plot_vectors", "plot_streamlines", "plot_field",
    "add_colorbar", "plot_contour",
    # data_extraction
    "profile_along_line", "tangents", "normals", "dist", "sort_indices",
    "sample_by_plane", "edge_lengths", "isoline",
    # quantities
    "momentum_thickness", "delta_star", "delta_99",
    # readers
    "Reader", "LegacyReader", "mark_boundary_cells", "NativeReader",
    "XMLReader", "VTUReader",
]

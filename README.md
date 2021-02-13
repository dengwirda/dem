## `RTopo-DEM`

Utilities for working with the `RTopo-2.x` data-sets and `MPAS-O`.

    python3 rtopo_map.py
    --mpas-file="path+name-to-mpas-mesh-file"
    --rtop-file="path+name-to-rtopo-pix-file"

can be used to inject (an unculled) `MPAS-O` mesh with the `bed_elevation`, `ocn_thickness` and `ice_thickness` varaibles derived from a compressed `RTopo-2.x` 'pixel' file (attached as assets to 'releases' of this repository). A 'remap'-style formulation is used, ensuring
consistent interpolation onto variable resolution voronoi-type grids.

    python3 rtopo_pix.py

can be used to create the compressed `RTopo-2.x` 'pixel' files used in the above workflow. The underlying `RTopo-2.x` data-sets are modified to use a cell-centred (rather than grid-centred) layout, compressed into `int16_t` arrays (to reduce memory), and are updated to expose the thickness-type variables listed above. To build the `RTopo-2.x` derived files, the full `60` and `30 arc-second` `RTopo-2.0.4` data-set must be downloaded: <a href="https://doi.pangaea.de/10.1594/PANGAEA.905295">`doi.pangaea.de/10.1594/PANGAEA.905295`</a>.

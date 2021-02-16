## `RTopo-DEM`

Utilities for working with the `RTopo-2.x`, `GEBCO` and `SRTM15+` data-sets in `MPAS-O`.

    python3 dem_map.py
    --mpas-file="path+name-to-mpas-mesh-file"
    --elev-file="path+name-to-DEM-pixel-file"

can be used to inject (an unculled) `MPAS-O` mesh with the `bed_elevation`, `ocn_thickness` and `ice_thickness` varaibles derived from a compressed `RTopo-2.x`, `GEBCO` or `SRTM15+` 'pixel' file (attached as assets to releases of this repository). A 'remap'-style formulation is used, ensuring consistent interpolation onto variable resolution voronoi-type grids.

    python3 dem_pix.py
    --elev-path="full-path-to-raw-DEM-assets"
    --save-path="full-path-to-output-storage"

can be used to create the compressed `RTopo-2.x`, `GEBCO` and `SRTM15+` 'pixel' files required for the above workflow. The underlying data-sets are modified to use a 'pixel'-centred (rather than grid-centred) layout, compressed into `int16_t` arrays (to reduce memory use), and are updated to expose the thickness-type variables listed above. In the case of `GEBCO` and `SRTM15+`, ice sheet/shelf elevation and thickness data is produced by 'blending' with the `R-Topo-2.x` data-sets at the ice fronts. To build the derived files, the following raw data-sets are needed: 

    * The `RTopo-2.0.4` data due to Schaffer et al: <a href="doi.pangaea.de/10.1594/PANGAEA.905295">`doi.pangaea.de/10.1594/PANGAEA.905295`</a>.
    * The `SRTM15+V2.1` data due to Tozer et al: <a href="doi.org/10.1029/2019EA000658">`doi.org/10.1029/2019EA000658`</a>.
    * The `GEBCO_v2020` data due to Weatherall et al: <a href="doi.org/10.5285/a29c5465-b138-234d-e053-6c86abc040b9">`doi.org/10.5285/a29c5465-b138-234d-e053-6c86abc040b9`</a>.

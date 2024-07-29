## `DEM`

Utilities for working with the `GEBCO` + `RTopo-2.x` digital elevation models in `MPAS`-like models.

    python3 dem_remap.py \
    --mpas-file="path+name-to-mpas-mesh-file" \
    --elev-file="path+name-to-DEM-pixel-file"

can be used to inject a (full-sphere) `MPAS`-like mesh with the `bed_elevation`, `ocn_thickness` and `ice_thickness` varaibles derived from a compressed `GEBCO` + `RTopo-2.x` 'pixel' file. A 'remap'-style formulation is used, ensuring consistent interpolation onto variable resolution voronoi-type grids.

    python3 dem_trnsf.py \
    --base-mesh="path+name-to-base-mesh-file" \
    --part-mesh="path+name-to-part-mesh-file"

can be used to transfer the remapped elevation data from a (full-sphere) `MPAS`-like mesh onto a partial sub-mesh, generated via a culling operation or equiv.

    python3 dem_pixel.py \
    --elev-path="full-path-to-raw-DEM-assets" \
    --save-path="full-path-to-output-storage"

can be used to create the compressed `GEBCO` + `RTopo-2.x` 'pixel' files required for the above workflow. The underlying data-sets are modified to use a 'pixel'-centred (rather than grid-centred) layout, compressed into `int16_t` arrays (to reduce memory use), and are updated to expose the thickness-type variables listed above. In the case of `GEBCO`, ice sheet/shelf elevation and thickness data is produced by 'blending' with the `R-Topo-2.x` data-sets at the ice fronts. To build the derived files, the following raw data-sets are needed: 

* The `RTopo-2.0.4` data due to Schaffer et al: <a href="doi.pangaea.de/10.1594/PANGAEA.905295">`doi.pangaea.de/10.1594/PANGAEA.905295`</a>.
* The `GEBCO_v2024` data due to Weatherall et al: <a href="doi.org/10.5285/1c44ce99-0a0d-5f4f-e063-7086abc0ea0f">`doi.org/10.5285/1c44ce99-0a0d-5f4f-e063-7086abc0ea0f`</a>.


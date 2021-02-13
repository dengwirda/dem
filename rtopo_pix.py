
import numpy as np
import netCDF4 as nc


def rtopo_60sec():
    """
    Create a zipped and pixel centred version of RTopo 2.0.4
    (60 arc-sec) to support remapping of elevation data.

    """
    # Authors: Darren Engwirda

    data = nc.Dataset("RTopo-2.0.4_1min_data.nc", "r")

    xpos = np.asarray(data["lon"][:], dtype=np.float64)
    ypos = np.asarray(data["lat"][:], dtype=np.float64)

    elev = np.asarray(
        data["bedrock_topography"][:], dtype=np.float32)
    surf = np.asarray(
        data["surface_elevation"][:], dtype=np.float32)
    base = np.asarray(
        data["ice_base_topography"][:], dtype=np.float32)

    elev = (elev[:-1:, :-1:] + elev[+1::, :-1:] + 
            elev[:-1:, +1::] + elev[+1::, +1::]) / 4.

    surf = (surf[:-1:, :-1:] + surf[+1::, :-1:] + 
            surf[:-1:, +1::] + surf[+1::, +1::]) / 4.

    base = (base[:-1:, :-1:] + base[+1::, :-1:] + 
            base[:-1:, +1::] + base[+1::, +1::]) / 4.

    elev = np.asarray(np.round(elev), dtype=np.int16)
    surf = np.asarray(np.round(surf), dtype=np.int16)
    base = np.asarray(np.round(base), dtype=np.int16)

    iceh = surf - base
    iceh[base == 0] = 0

    ocnh = np.maximum(0, base - elev)

    root = nc.Dataset(
        "RTopo_2_0_4_60sec_pixel.nc", "w", format="NETCDF4")
    root.description = "A zipped RTopo-2.0.4 (60 arc-sec) " + \
        "data-set, pixel centred and compressed via INT16."
    root.source = "RTopo-2.0.4_1min_data.nc"
    root.references = "doi.pangaea.de/10.1594/PANGAEA.905295"
    root.createDimension("num_lon", elev.shape[1] + 1)
    root.createDimension("num_col", elev.shape[1])    
    root.createDimension("num_lat", elev.shape[0] + 1)
    root.createDimension("num_row", elev.shape[0])

    data = root.createVariable("lon", "f8", ("num_lon"))
    data.units = "degrees_east"
    data[:] = xpos
    data = root.createVariable("lat", "f8", ("num_lat"))
    data.units = "degrees_north"
    data[:] = ypos
    data = root.createVariable(
        "bed_elevation", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = elev
    data = root.createVariable(
        "ice_thickness", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = iceh
    data = root.createVariable(
        "ocn_thickness", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = ocnh

    root.close()


def rtopo_30sec():
    """
    Create a zipped and pixel centred version of RTopo 2.0.4
    (30 arc-sec) to support remapping of elevation data.

    """
    # Authors: Darren Engwirda

    data = nc.Dataset(
        "RTopo-2.0.4_30sec_bedrock_topography.nc", "r")

    xpos = np.asarray(data["lon"][:], dtype=np.float64)
    ypos = np.asarray(data["lat"][:], dtype=np.float64)

    elev = np.asarray(
        data["bedrock_topography"][:], dtype=np.float32)

    data = nc.Dataset(
        "RTopo-2.0.4_30sec_surface_elevation.nc", "r")

    surf = np.asarray(
        data["surface_elevation"][:], dtype=np.float32)

    data = nc.Dataset(
        "RTopo-2.0.4_30sec_ice_base_topography.nc", "r")

    base = np.asarray(
        data["ice_base_topography"][:], dtype=np.float32)

    elev = (elev[:-1:, :-1:] + elev[+1::, :-1:] + 
            elev[:-1:, +1::] + elev[+1::, +1::]) / 4.

    surf = (surf[:-1:, :-1:] + surf[+1::, :-1:] + 
            surf[:-1:, +1::] + surf[+1::, +1::]) / 4.

    base = (base[:-1:, :-1:] + base[+1::, :-1:] + 
            base[:-1:, +1::] + base[+1::, +1::]) / 4.

    elev = np.asarray(np.round(elev), dtype=np.int16)
    surf = np.asarray(np.round(surf), dtype=np.int16)
    base = np.asarray(np.round(base), dtype=np.int16)

    iceh = surf - base
    iceh[base == 0] = 0

    ocnh = np.maximum(0, base - elev)

    root = nc.Dataset(
        "RTopo_2_0_4_30sec_pixel.nc", "w", format="NETCDF4")
    root.description = "A zipped RTopo-2.0.4 (30 arc-sec) " + \
        "data-set, pixel centred and compressed via INT16."
    root.source = "RTopo-2.0.4_30sec_data.nc"
    root.references = "doi.pangaea.de/10.1594/PANGAEA.905295"
    root.createDimension("num_lon", elev.shape[1] + 1)
    root.createDimension("num_col", elev.shape[1])    
    root.createDimension("num_lat", elev.shape[0] + 1)
    root.createDimension("num_row", elev.shape[0])

    data = root.createVariable("lon", "f8", ("num_lon"))
    data.units = "degrees_east"
    data[:] = xpos
    data = root.createVariable("lat", "f8", ("num_lat"))
    data.units = "degrees_north"
    data[:] = ypos
    data = root.createVariable(
        "bed_elevation", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = elev
    data = root.createVariable(
        "ice_thickness", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = iceh
    data = root.createVariable(
        "ocn_thickness", "i2", ("num_row", "num_col"))
    data.units = "m"
    data[:, :] = ocnh

    root.close()


if (__name__ == "__main__"):

    rtopo_60sec()
    rtopo_30sec()

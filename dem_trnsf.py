
import time
import numpy as np
import netCDF4 as nc
import argparse

# Authors: Darren Engwirda


def dem_trnsf(args):
    """
    Transfer elevation & ice+ocn-thickness data from a full 
    sphere MPAS mesh onto a partial sub-mesh resulting from
    a culling operation, or equiv.

    The sub-mesh is expected to be a piece of the base mesh
    such that (some) cells match 1:1.

    """

    base = nc.Dataset(args.base_mesh, "r")
    part = nc.Dataset(args.part_mesh, "r+")

#-- make a vector of cell-centre positions to match against

    xpos = np.vstack(( 
        np.vstack((
            base["xCell"][:],
            base["yCell"][:],
            base["zCell"][:])).T,
        np.vstack((
            part["xCell"][:],
            part["yCell"][:],
            part["zCell"][:])).T
    ))

    xidx = np.hstack((
        np.arange(
            0, base.dimensions["nCells"].size),
        np.arange(
            0, part.dimensions["nCells"].size)
    )) 

#-- culling shouldn't introduce fp round-off - but truncate 
#-- anyway...

    print("Build cell-to-cell maps...")
    
    ttic = time.time()

    xpos = np.round(xpos, decimals=8)

#-- use stable sorting to bring matching cell xyz (and idx)
#-- into "ascending" order

    imap = np.argsort(xpos[:, 2], kind="stable")
    xpos = xpos[imap, :]
    xidx = xidx[imap]
    imap = np.argsort(xpos[:, 1], kind="stable")
    xpos = xpos[imap, :]
    xidx = xidx[imap]
    imap = np.argsort(xpos[:, 0], kind="stable")
    xpos = xpos[imap, :]
    xidx = xidx[imap]

    diff = xpos[1:, :] - xpos[:-1, :]

    same = np.argwhere(np.logical_and.reduce((
        diff[:, 0] == 0, 
        diff[:, 1] == 0, 
        diff[:, 2] == 0))).ravel()

    ttoc = time.time()

    print("* built cell-to-cell map:",
          np.round(ttoc - ttic, decimals=1), "sec")

#-- cell inew in part matches iold in base - re-index elev.
#-- data-sets

    print("Transfer elevation data...")

    ttic = time.time()

    inew = xidx[same + 1]
    iold = xidx[same + 0]

#-- reindex cell elevation variables

    if ("bed_elevation" in base.variables.keys()):
        if ("bed_elevation" not in part.variables.keys()):
            part.createVariable("bed_elevation", "f4", ("nCells"))
    
        part["bed_elevation"][inew] = base["bed_elevation"][iold]

    if ("bed_slope" in base.variables.keys()):
        if ("bed_slope" not in part.variables.keys()):
            part.createVariable("bed_slope", "f4", ("nCells"))
    
        part["bed_slope"][inew] = base["bed_slope"][iold]

    if ("ocn_thickness" in base.variables.keys()):
        if ("ocn_thickness" not in part.variables.keys()):
            part.createVariable("ocn_thickness", "f4", ("nCells"))
    
        part["ocn_thickness"][inew] = base["ocn_thickness"][iold]

    if ("ice_thickness" in base.variables.keys()):
        if ("ice_thickness" not in part.variables.keys()):
            part.createVariable("ice_thickness", "f4", ("nCells"))
    
        part["ice_thickness"][inew] = base["ice_thickness"][iold]

#-- reindex cell histogram varaibles

    if ("nProfiles" in base.dimensions.keys()):
        if ("nProfiles" not in part.dimensions.keys()):
            part.createDimension(
                "nProfiles", base.dimensions["nProfiles"].size)

    if ("bed_elevation_profile" in base.variables.keys()):
        if ("bed_elevation_profile" not in part.variables.keys()):
            part.createVariable("bed_elevation_profile", 
                                "f4", ("nCells", "nProfiles"))

        part["bed_elevation_profile"][inew, :] = \
            base["bed_elevation_profile"][iold, :]

    if ("bed_slope_profile" in base.variables.keys()):
        if ("bed_slope_profile" not in part.variables.keys()):
            part.createVariable("bed_slope_profile", 
                                "f4", ("nCells", "nProfiles"))

        part["bed_slope_profile"][inew, :] = \
            base["bed_slope_profile"][iold, :]

    if ("ocn_thickness_profile" in base.variables.keys()):
        if ("ocn_thickness_profile" not in part.variables.keys()):
            part.createVariable("ocn_thickness_profile", 
                                "f4", ("nCells", "nProfiles"))

        part["ocn_thickness_profile"][inew, :] = \
            base["ocn_thickness_profile"][iold, :]

    if ("ice_thickness_profile" in base.variables.keys()):
        if ("ice_thickness_profile" not in part.variables.keys()):
            part.createVariable("ice_thickness_profile", 
                                "f4", ("nCells", "nProfiles"))

        part["ice_thickness_profile"][inew, :] = \
            base["ice_thickness_profile"][iold, :]
        
    ttoc = time.time()

    print("* reindex elevation data:",
          np.round(ttoc - ttic, decimals=1), "sec")

    base.close()
    part.close()

    return


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--base-mesh", dest="base_mesh", type=str,
        required=True, help="Name of (full) MPAS mesh.")

    parser.add_argument(
        "--part-mesh", dest="part_mesh", type=str,
        required=True, help="Name of culled MPAS mesh.")

    dem_trnsf(parser.parse_args())

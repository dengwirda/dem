
import os
import time
import numpy as np
import netCDF4 as nc
from scipy import spatial
from scipy import interpolate
from scipy.sparse import csr_matrix
import argparse


def map_to_r3(mesh, xlon, ylat, head, tail):
    """
    Map lon-lat coordinates to XYZ points. Restricted to the 
    panel LAT[HEAD:TAIL] to manage memory use.  

    """
    # Authors: Darren Engwirda

    sinx = np.sin(xlon * np.pi / 180.)
    cosx = np.cos(xlon * np.pi / 180.)
    siny = np.sin(ylat * np.pi / 180.)
    cosy = np.cos(ylat * np.pi / 180.)

    sinu, sinv = np.meshgrid(sinx, siny[head:tail])
    cosu, cosv = np.meshgrid(cosx, cosy[head:tail])

    rsph = mesh.sphere_radius

    xpos = rsph * cosu * cosv
    ypos = rsph * sinu * cosv
    zpos = rsph * sinv

    return np.vstack(
        (xpos.ravel(), ypos.ravel(), zpos.ravel())).T


def tria_area(rs, pa, pb, pc):
    """
    Calculate areas of spherical triangles [PA, PB, PC] on a
    sphere of radius RS.    

    """
    # Authors: Darren Engwirda

    lena = circ_dist(1., pa, pb)
    lenb = circ_dist(1., pb, pc)
    lenc = circ_dist(1., pc, pa)

    slen = 0.5 * (lena + lenb + lenc)

    tana = np.tan(0.5 * (slen - lena))
    tanb = np.tan(0.5 * (slen - lenb))
    tanc = np.tan(0.5 * (slen - lenc))

    edel = 4.0 * np.arctan(np.sqrt(
        np.tan(0.5 * slen) * tana * tanb * tanc))

    return edel * rs ** 2


def circ_dist(rs, pa, pb):
    """
    Calculate geodesic-length of great circles [PA, PB] on a
    sphere of radius RS.    

    """
    # Authors: Darren Engwirda    

    dlon = .5 * (pa[:, 0] - pb[:, 0])
    dlat = .5 * (pa[:, 1] - pb[:, 1])

    dist = 2. * rs * np.arcsin(np.sqrt(
        np.sin(dlat) ** 2 +
        np.sin(dlon) ** 2 * np.cos(pa[:, 1]) * np.cos(pb[:, 1])
    ))

    return dist


def cell_quad(mesh, ffun):
    """
    Eval. finite-volume integrals for mesh cells - splitting
    each into a series of triangles and employing an order-2
    quadrature rule.

    """
    # Authors: Darren Engwirda

    class base: pass

    abar = np.zeros(
        (mesh.dimensions["nCells"].size, 1), dtype=np.float64)
    fbar = np.zeros(
        (mesh.dimensions["nCells"].size, 1), dtype=np.float64)

    pcel = np.zeros(
        (mesh.dimensions["nCells"].size, 2), dtype=np.float64)
    pcel[:, 0] = mesh.variables["lonCell"][:]
    pcel[:, 1] = mesh.variables["latCell"][:]

    pcel = pcel * 180. / np.pi
    pcel[pcel[:, 0] > 180., 0] -= 360.

    fcel = ffun(pcel[:, 1], pcel[:, 0], grid=False)

    pcel = pcel * np.pi / 180.

    pvrt = np.zeros(
        (mesh.dimensions["nVertices"].size, 2), dtype=np.float64)
    pvrt[:, 0] = mesh.variables["lonVertex"][:]
    pvrt[:, 1] = mesh.variables["latVertex"][:]

    pvrt = pvrt * 180. / np.pi
    pvrt[pvrt[:, 0] > 180., 0] -= 360.

    fvrt = ffun(pvrt[:, 1], pvrt[:, 0], grid=False)

    pvrt = pvrt * np.pi / 180.

    cell = base()
    cell.edge = np.asarray(
        mesh.variables["edgesOnCell"][:], dtype=np.int32)
    cell.topo = np.asarray(
        mesh.variables["nEdgesOnCell"][:], dtype=np.int32)

    edge = base()
    edge.vert = np.asarray(
        mesh.variables["verticesOnEdge"][:], dtype=np.int32)
    
    for epos in range(np.max(cell.topo)):

        mask = cell.topo > epos

        icel = np.argwhere(mask).ravel()

        ifac = cell.edge[mask, epos] - 1

        ivrt = edge.vert[ifac, 0] - 1
        jvrt = edge.vert[ifac, 1] - 1

        rsph = mesh.sphere_radius

        atri = tria_area(
            rsph, pcel[icel], pvrt[ivrt], pvrt[jvrt])

        atri = np.reshape(atri, (atri.size, 1))

        ftri = (fcel[icel] + fvrt[ivrt] + fvrt[jvrt])

        ftri = np.reshape(ftri, (ftri.size, 1))

        abar[icel] += atri
        fbar[icel] += atri * ftri / 3.

    return fbar / abar


def rtopo_map(args):
    """
    Map elevation and ice+ocn-thickness data from a "zipped" 
    RTopo data-set onto the cells in an MPAS mesh.

    Cell values are a blending of an approx. integral remap
    and a local quadrature rule.

    """
    # Authors: Darren Engwirda

    print("Loading assests...")

    rtop = nc.Dataset(args.rtop_file, "r")
    mesh = nc.Dataset(args.mpas_file, "r+")

#-- Compute an approximate remapping, associating pixels in
#-- the DEM with cells in the MPAS mesh. Since polygons are
#-- Voronoi, the point-in-cell query can be computed by
#-- finding nearest neighbours. This remapping is an approx.
#-- as no partial pixel-cell intersection is computed.  

    print("Building KDtree...")
 
    ppos = np.zeros(
        (mesh.dimensions["nCells"].size, 3), dtype=np.float64)
    ppos[:, 0] = mesh["xCell"][:]
    ppos[:, 1] = mesh["yCell"][:]
    ppos[:, 2] = mesh["zCell"][:]

    tree = spatial.cKDTree(ppos, leafsize=8)

    print("Remap elevation...")

    xlon = np.asarray(rtop["lon"][:], dtype=np.float64)
    ylat = np.asarray(rtop["lat"][:], dtype=np.float64)

    xmid = .5 * (xlon[:-1:] + xlon[1::])
    ymid = .5 * (ylat[:-1:] + ylat[1::])

    indx = np.asarray(np.round(
        np.linspace(-1, ymid.size, 9)), dtype=np.int32)

    print("* process over 8 tiles.")

    nset = []
    for tile in range(indx.size - 1):

        head = indx[tile + 0] + 1
        tail = indx[tile + 1] + 1

        qpos = map_to_r3(mesh, xmid, ymid, head, tail)

        ttic = time.time()
        __, nloc = tree.query(qpos, n_jobs=-1)
        ttoc = time.time()
        print("* built node-to-cell map:", ttoc - ttic)
        
        nset.append(nloc)

    near = np.concatenate(nset)
    
    del tree; del qpos; del nset; del nloc

#-- Build cell-to-pixel map as a sparse matrix, and average
#-- RTopo pixel values within each cell.
    
    cols = np.arange(0, near.size)
    vals = np.ones(near.size, dtype=np.int8)

    smat = csr_matrix((vals, (near, cols)))
    
    nmap = np.sum(smat, axis=1)
    
    vals = np.asarray(
        rtop["bed_elevation"][:], dtype=np.float32)
    vals = np.reshape(vals, (vals.size, 1))

    emap = (smat * vals) / np.maximum(+1, nmap)

    vals = np.asarray(
        rtop["ocn_thickness"][:], dtype=np.float32)
    vals = np.reshape(vals, (vals.size, 1))

    omap = (smat * vals) / np.maximum(+1, nmap)

    vals = np.asarray(
        rtop["ice_thickness"][:], dtype=np.float32)
    vals = np.reshape(vals, (vals.size, 1))

    imap = (smat * vals) / np.maximum(+1, nmap)

    del smat; del cols; del vals; del near

#-- If the resolution of the mesh is greater, or comparable 
#-- to that of the DEM, the approx. remapping (above) will 
#-- result in a low order interpolation.
#-- Thus, blend with a local order-2 quadrature formulation

    print("Eval. elevation...")

    vals = np.asarray(
        rtop["bed_elevation"][:], dtype=np.float32)

    ttic = time.time()
    eint = cell_quad(mesh, 
        interpolate.RectBivariateSpline(ymid, xmid, vals, kx=1, ky=1))
    ttoc = time.time()
    print("* compute cell integrals:", ttoc - ttic)

    vals = np.asarray(
        rtop["ocn_thickness"][:], dtype=np.float32)

    ttic = time.time()
    oint = cell_quad(mesh, 
        interpolate.RectBivariateSpline(ymid, xmid, vals, kx=1, ky=1))
    ttoc = time.time()
    print("* compute cell integrals:", ttoc - ttic)

    vals = np.asarray(
        rtop["ice_thickness"][:], dtype=np.float32)

    ttic = time.time()
    iint = cell_quad(mesh, 
        interpolate.RectBivariateSpline(ymid, xmid, vals, kx=1, ky=1))
    ttoc = time.time()
    print("* compute cell integrals:", ttoc - ttic)
    

    print("Save to dataset...")
    
    ebar = (np.multiply(nmap, emap) + eint) / (1 + nmap)
    obar = (np.multiply(nmap, omap) + oint) / (1 + nmap)
    ibar = (np.multiply(nmap, imap) + iint) / (1 + nmap)
    
    if ("bed_elevation" in mesh.variables.keys()):
        mesh["bed_elevation"][:] = ebar
    else:
        data = mesh.createVariable("bed_elevation", "f8", ("nCells"))
        data[:] = ebar

    if ("ocn_thickness" in mesh.variables.keys()):
        mesh["ocn_thickness"][:] = obar
    else:
        data = mesh.createVariable("ocn_thickness", "f8", ("nCells"))
        data[:] = obar

    if ("ice_thickness" in mesh.variables.keys()):
        mesh["ice_thickness"][:] = ibar
    else:
        data = mesh.createVariable("ice_thickness", "f8", ("nCells"))
        data[:] = ibar

    mesh.close()


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--mpas-file", dest="mpas_file", type=str,
        required=True, help="Path to user MPAS file.")

    parser.add_argument(
        "--rtop-file", dest="rtop_file", type=str,
        required=True, help="Path to RTopo PIX file.")

    rtopo_map(parser.parse_args())

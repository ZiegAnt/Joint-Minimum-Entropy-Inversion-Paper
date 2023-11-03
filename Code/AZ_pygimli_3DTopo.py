# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:11:09 2023

@author: azieg
"""
# %% IMPORTS
import numpy             as np
import pygimli           as pg
import pygimli.meshtools as mt
import functools

def AZ_createParaMeshSurface(sensors,boundary=[0,0],surfaceMeshQuality=30, addTopo=None, boundaryMarker=None):
    '''

    Parameters
    ----------
    sensors : np.asarray([x, y, z]).T
        points that cover area
    boundary : list
        that is appended to sensor area. The default is [0,0].
    surfaceMeshQuality : int, optional
        The default is 30.
    surfaceMeshArea : int, optional
        The default is 0.
    addTopo : list, optional
        3D topography points. The default is None.

    Returns
    -------
    s : surface mesh
        3D surface mesh.

    '''

    if hasattr(sensors, 'sensors'):
        sensors = sensors.sensors()
    sensors = np.asarray(sensors)


    if boundary is None:
        boundary = [0.0, 0.0]

    # find maximum extent
    maxx = max(sensors[:,0])
    maxy = max(sensors[:,1])
    minx = min(sensors[:,0])
    miny = min(sensors[:,1])
    boundaryRect = mt.createRectangle(start=[minx-boundary[0],maxy+boundary[1]],
                                                end=[maxx+boundary[0],miny-boundary[1]])
    
    for i in range(4):
        boundaryRect.boundary(i).setMarker(i+1)
        boundaryRect.node(i).setMarker((i % 4 + 1) * 10)

    # collect all pnts with topography
    if addTopo is not None:
        if min(sensors[:, 2]) != max(sensors[:, 2]) and sensors[0][2] != 0.0:
            pnts = np.vstack((sensors, addTopo))
        else:
            pnts = np.asarray(addTopo)
    else:
        pnts = np.array(sensors)

    # add maximal extent corners to median topo
    boundaryRect.translate([0, 0, np.median(pnts[:, 2])])
    pnts = np.vstack((boundaryRect.positions(), pnts))

    # create mesh for topo interpolation
    pntsSurface = mt.createMesh(pnts[:, 0:2])

    # find parameter extent
    paraRect = mt.createRectangle(start=[minx,maxy],
                                            end=[maxx,miny])


    # create surface mesh of sensors and with maximal and parameter extent
    surfacePLC = boundaryRect + pnts #+ sensors[:, 0:2] 
    surface = mt.createMesh(surfacePLC, quality=surfaceMeshQuality)


    # interpolate Topography to surface
    sZ = pg.interpolate(pntsSurface, pnts[:, 2], surface.positions())
    for n in surface.nodes():
        n.translate(0, 0, sZ[n.id()])

#     # create 3D surfacemesh
#     if surface.dimension() != 2:
#         pg.error("Need two dimensional mesh")
#     if surface.cellCount() == 0:
#         pg.error("Need a two dimensional mesh with cells")

#     s = pg.Mesh(dim=3, isGeometry=True)

#     idx=[s.createNode(n.pos(), n.marker()).id() for n in surface.nodes()]

#     for c in surface.cells():
#         if c.ids()[0]<len(s.nodes()) and c.ids()[1]<len(s.nodes()) and c.ids()[2]<len(s.nodes()):
#             s.createBoundary(c.ids(), marker=c.marker())

#     if boundaryMarker is not None:
#         s.setBoundaryMarkers(np.full(s.boundaryCount(),
#                                             boundaryMarker))

     # create 3D surfacemesh
    s = pg.meshtools.createSurface(surface, boundaryMarker=pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
        
    return s



def AZ_createParaMeshPLC3D(sensors, paraDX=0, paraDepth=-1, boundary=None,surfaceMeshQuality=30, addTopo=None, **kwargs):
    """Create a geometry (PLC) for an 3D inversion parameter mesh.


        Args
        ----
        sensors: Sensor list or pg.DataContainer with .sensors()
            Sensor positions.

        paraDX : float [1]
            Absolute distance for node refinement (0=none).
            Refinement node will be placed below the surface.

        paraDepth : float[-1], optional
            Maximum depth in m for parametric domain.
            Automatic (<=0) results in 0.4 * maximum sensor span range in m.
            Depth is set to median sensors depth + paraDepth.

        boundary: [float, float] [0., 0.]
            Boundary width to be appended for domain prolongation in relative
            para domain size.

        surfaceMeshQuality: float [30]
            Quality of the surface mesh.

        addTopo: [[x,y,z],]
            Number of additional nodes for topography.

        Returns
        -------
        poly: :gimliapi:`GIMLI::Mesh`
            Piecewise linear complex (PLC) containing nodes and edges
    """

    
    if hasattr(sensors, 'sensors'):
        sensors = sensors.sensors()

    sensors = np.asarray(sensors)

    if boundary is None:
        boundary = [0.0, 0.0]

    surface = AZ_createParaMeshSurface(sensors,
                                       boundary=boundary,
                                       surfaceMeshQuality=surfaceMeshQuality,
                                       addTopo=addTopo, boundaryMarker=-2)

    # find depth and paradepth
    xSpan = (max(sensors[:, 0]) - min(sensors[:, 0]))
    ySpan = (max(sensors[:, 1]) - min(sensors[:, 1]))

    if paraDepth == -1:
        paraDepth = (0.4*(max(xSpan, ySpan)))

    paraDepth = np.median(sensors[:, 2]) - paraDepth
    depth = paraDepth #- max(boundary[0]*xSpan, boundary[1]*ySpan)/2

    def sortP(p):
        base = pg.core.Line(p[0], p[1]).at(-1e7)

        def cmp_(p1, p2):
            if p1.distSquared(base) < p2.distSquared(base):
                return -1
            else:
                return 1

        p.sort(key=functools.cmp_to_key(cmp_))

    bounds = [surface]
    # close outer surfaces
    bttm = []
    for i in range(4):

        p = [n.pos() for n in surface.nodes() if n.marker() == i+1]
        p.append(surface.nodes(
            surface.nodeMarkers() == (i % 4 + 1) * 10)[0].pos())
        p.append(surface.nodes(
            surface.nodeMarkers() == ((i + 1) % 4 + 1) * 10)[0].pos())
        sortP(p)

        p0 = pg.Pos(p[-1])
        p0[2] = depth
        p.append(p0)

        p1 = pg.Pos(p[0])
        p1[2] = depth
        p.append(p1)

        m = mt.createPolygon(p, isClosed=True)
        # f = mt.createFacet(m, boundaryMarker=-2)
        f = mt.createFacet(m, boundaryMarker=0)

        bounds.append(f)
        bttm.append(p0)
        bttm.append(p1)

    m = mt.createRectangle(pnts=bttm, minBB=True)
    m.translate([0, 0, depth])

    bttmA = mt.createFacet(m, boundaryMarker=1)
    bounds.append(bttmA)

    pdPLC = mt.mergePLC(bounds)

    if paraDX > 0:
        for s in sensors:
            pdPLC.createNode(s - [0.0, 0.0, paraDX])

    pdPLC.addRegionMarker(pg.center(bttmA.positions()) + [0.0, 0.0, 0.1],
                          marker=1)

    return pdPLC
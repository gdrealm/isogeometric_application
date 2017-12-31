import math
#from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

###
### This module is a factory to generate typical geometries for isogeometric analysis, e.g. circle, l-shape, ...
###

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
bsplines_patch_util = BSplinesPatchUtility()

### Create a line from start_point to end_point with knot vector [0 0 0 ... 1 1 1]
### On output the pointer to the patch will be returned
def CreateLine(start_point, end_point, order = 1):
    Id = 0
    fes = nurbs_fespace_library.CreateLinearFESpace(order)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(start_point[0], start_point[1], start_point[2], fes.Number(0), end_point[0], end_point[1], end_point[2])
    patch_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create an arc at center on the surface perpendicular to the axis. By default, the quadratic arc is generated. The knot vector will be [0 0 0 1 1 1]
### On output the pointer to the patch will be returned. Small arc means that the open angle is less than 90 degrees.
def CreateSmallArc(center, axis, radius, start_angle, end_angle):
    ## firstly create an arc in xy plane at (0, 0)
    Id = 0
    fes = nurbs_fespace_library.CreateLinearFESpace(2)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes.Number(0), radius, 0.0, 0.0)

    sweep = end_angle - start_angle
    if abs(sweep > 90):
        raise ValueError('the open angle must be in [-90, 90] degrees, sweep =', sweep)

    dsweep = 0.5*sweep/180.0*math.pi
    wm = math.cos(dsweep)
    x = radius*wm
    y = radius*math.sin(dsweep)
    xm = x + y*math.tan(dsweep)

    if axis == 'z':
        trans = RotationZ(start_angle + 0.5*sweep)
    elif axis == 'y':
        trans = RotationZ(start_angle + 0.5*sweep)
        trans.AppendTransformation(RotationX(90.0))
    elif axis == 'x':
        trans = RotationZ(start_angle + 0.5*sweep + 90.0)
        trans.AppendTransformation(RotationY(90.0))
    trans.AppendTransformation(Translation(center[0], center[1], center[2]))

    pt1 = ctrl_grid[0]
    pt1.WX = x
    pt1.WY = -y
    pt1.WZ = 0.0
    pt1.W = 1.0
    pt1.ApplyTransformation(trans)
    ctrl_grid[0] = pt1

    pt2 = ctrl_grid[1]
    pt2.WX = wm*xm
    pt2.WY = 0.0
    pt2.WZ = 0.0
    pt2.W = wm
    pt2.ApplyTransformation(trans)
    ctrl_grid[1] = pt2

    pt3 = ctrl_grid[2]
    pt3.WX = x
    pt3.WY = y
    pt3.WZ = 0.0
    pt3.W = 1.0
    pt3.ApplyTransformation(trans)
    ctrl_grid[2] = pt3

    patch_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a ring at center on the surface perpendicular to the axis. By default, the quadratic arc is generated. The knot vector will be [0 0 0 1 1 1]
### On output the pointer to the patch will be returned. Small ring means that the open angle is less than 90 degrees.
def CreateSmallRing(center, axis, rin, rout, start_angle, end_angle):
    ## create inner arc
    iarc_ptr = CreateSmallArc(center, axis, rin, start_angle, end_angle)
    iarc = iarc_ptr.GetReference()

    ## create outer arc
    oarc_ptr = CreateSmallArc(center, axis, rout, start_angle, end_angle)
    oarc = oarc_ptr.GetReference()

    ## create ring
    ring_patch_ptr = bsplines_patch_util.CreateConnectedPatch(iarc, oarc)
    return ring_patch_ptr


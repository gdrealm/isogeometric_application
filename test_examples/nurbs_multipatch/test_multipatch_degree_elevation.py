##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_util = ControlGridUtility()
multipatch_util = MultiPatchUtility()
mpatch = MultiPatch2D()

def CreateMultiPatch():

    fes1 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_1 = grid_util.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
    patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
    patch1 = patch1_ptr.GetReference()
    patch1.CreateControlPointGridFunction(ctrl_grid_1)
    #print(patch1)
    multipatch_util.ExportGeo(patch1, "patch1.txt")

    fes2 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_2 = grid_util.CreateRectangularControlPointGrid(1.0, 0.0, fes1.Number(0), fes1.Number(1), 2.0, 1.0)
    patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
    #print(patch2)
    multipatch_util.ExportGeo(patch2, "patch2.txt")

    mpatch.AddPatch(patch1)
    mpatch.AddPatch(patch2)
    mpatch.ResetId()
    mpatch.MakeNeighbor(patch1, BoundarySide.Right, patch2, BoundarySide.Left)

    print("############REFINEMENT###############")
    multipatch_refine_util = MultiPatchRefinementUtility()

    multipatch_refine_util.DegreeElevate(patch1_ptr, [0, 1])

    #patch1 = patch1_ptr.GetReference()
    #patch2 = patch2_ptr.GetReference()

    return mpatch


#print("############RESULTS###############")
mpatch = CreateMultiPatch()
print(mpatch)
multipatch_util.ExportGlvis(mpatch, "mpatch.mesh")


#################RESULTS#####################(Validated with Matlab)
##>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<
##-------------Begin MultiPatchInfo-------------
##MultiPatch overview: Number of patches = 2
##MultiPatch details:
##-------------Begin PatchInfo-------------
##Patch2D, Id = 1, Add = 0x12d56d0
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x12d4f90, n = ( 5 5), p = ( 4 4)
## knot vector 0: 0 0 0 0 0 1 1 1 1 1
## knot vector 1: 0 0 0 0 0 1 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [5, 5]
## Data:
## (
##  ( (0, 0, 0, 1) (0, 0.25, 0, 1) (0, 0.5, 0, 1) (0, 0.75, 0, 1) (0, 1, 0, 1))
##  ( (0.25, 0, 0, 1) (0.25, 0.25, 0, 1) (0.25, 0.5, 0, 1) (0.25, 0.75, 0, 1) (0.25, 1, 0, 1))
##  ( (0.5, 0, 0, 1) (0.5, 0.25, 0, 1) (0.5, 0.5, 0, 1) (0.5, 0.75, 0, 1) (0.5, 1, 0, 1))
##  ( (0.75, 0, 0, 1) (0.75, 0.25, 0, 1) (0.75, 0.5, 0, 1) (0.75, 0.75, 0, 1) (0.75, 1, 0, 1))
##  ( (1, 0, 0, 1) (1, 0.25, 0, 1) (1, 0.5, 0, 1) (1, 0.75, 0, 1) (1, 1, 0, 1))
## )
##Neighbors =  right:2
##-------------End PatchInfo-------------
##-------------Begin PatchInfo-------------
##Patch2D, Id = 2, Add = 0x12d4e20
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x135df50, n = ( 4 5), p = ( 3 4)
## knot vector 0: 0 0 0 0 1 1 1 1
## knot vector 1: 0 0 0 0 0 1 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [4, 5]
## Data:
## (
##  ( (1, 0, 0, 1) (1, 0.25, 0, 1) (1, 0.5, 0, 1) (1, 0.75, 0, 1) (1, 1, 0, 1))
##  ( (1.33333, 0, 0, 1) (1.33333, 0.25, 0, 1) (1.33333, 0.5, 0, 1) (1.33333, 0.75, 0, 1) (1.33333, 1, 0, 1))
##  ( (1.66667, 0, 0, 1) (1.66667, 0.25, 0, 1) (1.66667, 0.5, 0, 1) (1.66667, 0.75, 0, 1) (1.66667, 1, 0, 1))
##  ( (2, 0, 0, 1) (2, 0.25, 0, 1) (2, 0.5, 0, 1) (2, 0.75, 0, 1) (2, 1, 0, 1))
## )
##Neighbors =  left:1
##-------------End PatchInfo-------------
##-------------End MultiPatchInfo-------------
##>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<



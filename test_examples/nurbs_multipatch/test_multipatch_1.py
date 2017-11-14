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

nurbs_patch_library = NURBSPatchLibrary()
grid_func_util = GridFunctionUtility()
nurbs_patch_util = NURBSPatchUtility()
mpatch = MultiPatch2D()

patch1_ptr = nurbs_patch_library.CreateRectangularPatch(3, 3)
patch1 = patch1_ptr.GetReference()
#print("patch1 address:")
#mpatch.PrintAddress(patch1)

#patch1 = nurbs_patch_library.CreateRectangularPatch(3, 3)
ctrl_grid_1 = grid_func_util.CreateRectangularControlPointGrid(0.0, 0.0, patch1.Number(0), patch1.Number(1), 1.0, 1.0)
patch1.SetControlPointGrid(ctrl_grid_1)
print(patch1)
nurbs_patch_util.ExportGeo(patch1, "patch1.txt")

patch2_ptr = nurbs_patch_library.CreateRectangularPatch(3, 3)
patch2 = patch2_ptr.GetReference()
#print("patch2 address:")
#mpatch.PrintAddress(patch2)
#patch2 = nurbs_patch_library.CreateRectangularPatch(3, 3)
ctrl_grid_2 = grid_func_util.CreateRectangularControlPointGrid(1.0, 1.0, patch2.Number(0), patch2.Number(1), 2.0, 2.0)
patch2.SetControlPointGrid(ctrl_grid_2)
print(patch2)
nurbs_patch_util.ExportGeo(patch2, "patch2.txt")

mpatch.AddPatch(patch1)
mpatch.AddPatch(patch2)
mpatch.ResetId()
mpatch.MakeNeighbor(patch1, BoundarySide.Right, patch2, BoundarySide.Left)
#print(mpatch)

print("############REFINEMENT###############")
multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.InsertKnots(patch1_ptr, [[0.5], [0.5]])
#patch1 = patch1_ptr.GetReference()
#patch2 = patch2_ptr.GetReference()
#print("############RESULTS###############")
print(mpatch)


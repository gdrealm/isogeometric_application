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

patch1_ptr = nurbs_patch_library.CreateLinearPatch(3)
patch1 = patch1_ptr.GetReference()

ctrl_grid_1 = grid_func_util.CreateLinearControlPointGrid(0.0, 0.0, 0.0, patch1.Number(0), 1.0, 0.0, 0.0)
patch1.SetControlPointGrid(ctrl_grid_1)

print("######BEFORE REFINEMENT#######")
print(patch1)
nurbs_patch_util = NURBSPatchUtility()
nurbs_patch_util.ExportGeo(patch1, "line.txt")

mpatch = MultiPatch1D()
mpatch.AddPatch(patch1)
mpatch.ResetId()

multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.DegreeElevate(patch1_ptr, [1])
print("######AFTER REFINEMENT#######")
print(mpatch) # result is the same as matlab code below

##lin = nrbdegelev(nrbline,2)
##lin1 = nrbdegelev(lin,1)
##lin1.coefs


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

patch1 = nurbs_patch_library.CreateRectangularPatch(3, 3)
ctrl_grid_1 = grid_func_util.CreateRectangularControlPointGrid(0.0, 0.0, patch1.Number(0), patch1.Number(1), 1.0, 1.0)
patch1.SetControlPointGrid(ctrl_grid_1)
print(patch1)

patch2 = nurbs_patch_library.CreateCubicPatch(3, 3, 2)
ctrl_grid_2 = grid_func_util.CreateCubicControlPointGrid(0.0, 0.0, 0.0, patch2.Number(0), patch2.Number(1), patch2.Number(2), 1.0, 1.0, 1.0)
patch2.SetControlPointGrid(ctrl_grid_2)
print(patch2)


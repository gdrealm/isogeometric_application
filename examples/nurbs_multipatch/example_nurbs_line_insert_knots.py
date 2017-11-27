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
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()

fes1 = nurbs_fespace_library.CreateLinearFESpace(3)
ctrl_grid_1 = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes1.Number(0), 1.0, 0.0, 0.0)
#print(fes1)
#print(ctrl_grid_1)

patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
print(patch1)

print("######BEFORE REFINEMENT#######")
print(patch1)
multipatch_util.ExportGeo(patch1, "line.txt")

print("######AFTER REFINEMENT#######")
mpatch = MultiPatch1D()
mpatch.AddPatch(patch1)
mpatch.ResetId()

multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.InsertKnots(patch1_ptr, [[0.4, 0.5]])
print(mpatch) # result is the same as matlab code below
multipatch_util.ExportGlvis(mpatch, "line.mesh")

##lin = nrbdegelev(nrbline,2)
##lin1 = nrbkntins(lin,[0.4 0.5])
##lin1.coefs



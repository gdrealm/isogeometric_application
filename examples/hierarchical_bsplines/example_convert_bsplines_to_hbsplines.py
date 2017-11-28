##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()
hbsplines_patch_util = HBSplinesPatchUtility()

fn = "../../tests/infinite_plate.txt"
patches_ptr = bsplines_patch_util.CreatePatchFromGeo(fn)
patch = patches_ptr[0].GetReference()
print(patch)

hpatch_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch)
hpatch = hpatch_ptr.GetReference()
print(hpatch)



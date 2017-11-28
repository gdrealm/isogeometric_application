##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()
hbsplines_patch_util = HBSplinesPatchUtility()
hbsplines_refinement_util = HBSplinesRefinementUtility()

def CreatePatch():
    fn = "square.txt"
    patches_ptr = bsplines_patch_util.CreatePatchFromGeo(fn)
    patch = patches_ptr[0].GetReference()
    #print(patch)

    hpatch_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch)
    hpatch = hpatch_ptr.GetReference()
    #print(hpatch)

    echo_level = 1
    hbsplines_refinement_util.Refine(hpatch, 1, echo_level)
    hpatch.FESpace().Enumerate()
#    print(hpatch)

    return hpatch_ptr

def main():
    patch_ptr = CreatePatch()
    patch = patch_ptr.GetReference()
    print(patch)

if __name__ == "__main__":
    main()


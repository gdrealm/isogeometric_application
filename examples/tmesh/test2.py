##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

## test the tsplines bezier extraction functionality
# utils = BezierUtils()
# utils.test_tsplines_1() # test the compute_extended_knot_vector


## test the functionality of the Tmesh
Tmesh = TsMesh2D()

Tmesh.BeginConstruct()
Tmesh.ReadFromFile('test2.tmesh')
Tmesh.EndConstruct()
# print(Tmesh)

Tmesh.BuildExtendedTmesh()
print "Tmesh is analysis suitable: ", Tmesh.IsAnalysisSuitable()

# Tmesh.ExportMatlab("tmesh2.m", "topology")
# Tmesh.ExportMatlab("tmesh2.m", "knots")

Tmesh.BuildAnchors('test2.coordinates')

Tmesh.BuildCells()

Tmesh.ExportMDPA('tmesh2.mdpa', 1, 1)


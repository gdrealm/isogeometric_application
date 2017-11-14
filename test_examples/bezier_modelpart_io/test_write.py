##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

HNmesh = HnMesh("My Mesh")
HNmesh.ReadMesh("geo_ring.txt")
HNmesh.SetEchoLevel(HN_ECHO_FLAGS.ECHO_REFIMENT)
HNmesh.Refine(1)
HNmesh.LinearDependencyRefine(1)
#HNmesh.CheckNestedSpace()
HNmesh.PrintRefinementHistory()

print(HNmesh)
HNmesh.BuildMesh()
HNmesh.PrintKnotVectors()
# HNmesh.ExportCellTopology("/home/hbui/workspace2/isogeometric/hierarchical_b-splines/geo_ring_hnurbs_cell_topo.m", False)
# HNmesh.ExportMatlab("/home/hbui/workspace2/isogeometric/hierarchical_b-splines/geo_ring_hnurbs.m")
# HNmesh.ExportMDPA("geo_ring_hnurbs.mdpa")
HNmesh.ExportMDPA2("geo_ring_hnurbs.mdpa")
# HNmesh.ExportPostMDPA("geo_ring_hnurbs_post.mdpa", 31, 31, 11)


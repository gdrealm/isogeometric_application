##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

model_part_io = BezierModelPartIO("model")
mp = ModelPart()

model_part_io.ReadModelPart(mp)


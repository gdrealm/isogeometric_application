##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

dm = DomainManager2D(1)
dm.AddXcoord(0.0)
dm.AddXcoord(0.5)
dm.AddXcoord(1.0)
dm.AddXcoord(1.5)
dm.AddXcoord(2.0)
dm.AddYcoord(0.0)
dm.AddYcoord(0.5)
dm.AddYcoord(1.0)
dm.AddYcoord(1.5)
dm.AddYcoord(2.0)
dm.AddYcoord(3.0)

dm.AddCell(0.0, 0.5, 0.0, 0.5)
dm.AddCell(0.0, 0.5, 0.5, 1.0)
dm.AddCell(0.5, 1.0, 0.0, 0.5)
dm.AddCell(0.5, 1.0, 0.5, 1.0)
dm.AddCell(1.5, 2.0, 1.5, 2.0)

print(dm)

# is_inside = dm.IsInside(0.25, 0.75, 0.25, 0.75) # True, correct
# is_inside = dm.IsInside(0.25, 0.5, 0.25, 0.5) # True, correct
# is_inside = dm.IsInside(0.25, 0.5, 0.25, 1.1) # False, correct
# is_inside = dm.IsInside(-0.25, 0.5, 0.25, 1.0) # False, correct
# is_inside = dm.IsInside(-0.0, 0.5, 0.25, 1.0) # True, correct
# is_inside = dm.IsInside(0.0, 1.0, 0.0, 1.0) # True, correct
# is_inside = dm.IsInside(-0.001, 1.0, 0.0, 1.0) # False, correct
# is_inside = dm.IsInside(0.8, 1.6, 0.8, 1.6) # False, correct
is_inside = dm.IsInside(1.6, 1.8, 1.6, 1.8) # True, correct
print(is_inside)





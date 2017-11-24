##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your actual configuration                #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./geo_ring')
import geo_ring_include
from geo_ring_include import *

# Initialize model
model = geo_ring_include.Model('geo_ring', os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF infinite_plate.gid |
# vvvvvvvvvvvvvvvvvvvvv

time = 0.0
model.Solve(time, 0, 0, 0, 0)
model.WriteOutput(time)

#print model.model_part.Elements[1]

test_utils = NURBSTestUtils()
#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 0.0, 0.0)
#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 0.0, 1.0)
#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 1.0, 0.0)
#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 1.0, 1.0)

#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 0.5, 0.5)
#test_utils.ProbeShapeFunctionValues(model.model_part.Elements[1], 0.5, 0.5)
#test_utils.ProbeShapeFunctionDerivatives(model.model_part.Elements[1], 0.5, 0.5)
#test_utils.ProbeJacobian(model.model_part.Elements[1], 0.5, 0.5)

#test_utils.ProbeGlobalCoordinates(model.model_part.Elements[1], 0.0, 0.5)
#test_utils.ProbeShapeFunctionValues(model.model_part.Elements[1], 0.0, 0.5)
#test_utils.ProbeShapeFunctionDerivatives(model.model_part.Elements[1], 0.0, 0.5)
#test_utils.ProbeJacobian(model.model_part.Elements[1], 0.0, 0.5)
#test_utils.ProbeJacobian(model.model_part.Elements[1], 0.0, 0.0)

test_utils.ProbeShapeFunctionValues(model.model_part.Elements[1], 0.333333, 0.333333)


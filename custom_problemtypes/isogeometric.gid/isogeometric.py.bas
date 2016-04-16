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
sys.path.append('./rEpLaCeMeNtStRiNg.gid')
import rEpLaCeMeNtStRiNg_include
from rEpLaCeMeNtStRiNg_include import **

# Initialize model
model = rEpLaCeMeNtStRiNg_include.Model('rEpLaCeMeNtStRiNg', os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF rEpLaCeMeNtStRiNg.gid |
# vvvvvvvvvvvvvvvvvvvvv


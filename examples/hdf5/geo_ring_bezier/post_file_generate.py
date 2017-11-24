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
###  POST-PROCESSING  ############################################
##################################################################

## create model_part post # this is done for only one time
model_part_post = ModelPart("isogeometric_mesh")
model_part_post.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part_post.AddNodalSolutionStepVariable(REACTION)
model_part_post.AddNodalSolutionStepVariable(STRESSES)
import structural_solver_advanced
structural_solver_advanced.AddDofs(model_part_post)
#model.isogeometric_post_utility.GenerateModelPart(model_part_post, PostElementType.Quadrilateral)
#model.isogeometric_post_utility.GenerateModelPart2(model_part_post)
model.isogeometric_post_utility.GenerateModelPart2AutoCollapse(model_part_post, 0.3, 0.3, 0.3, 1.0e-3)
########################################

## create the gid_io
write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
write_elements = WriteConditionsFlag.WriteConditions
#write_elements = WriteConditionsFlag.WriteElementsOnly
post_mode = GiDPostMode.GiD_PostAscii
multi_file_flag = MultiFileFlag.MultipleFiles
gid_io = StructuralGidIO(model.path+model.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements)
########################################


## read data
time = 0.0
model.model_part.CloneTimeStep(time)
hdf5_post = HDF5PostUtility(model.problem_name + "_" + str(time) + ".h5", "Read-Only")
hdf5_post.ReadNodalResults(DISPLACEMENT, model.model_part)
hdf5_post.ReadNodalResults(REACTION, model.model_part)
hdf5_post.ReadNodalResults(STRESSES, model.model_part)

for node in model.model_part.Nodes:
    print node.GetSolutionStepValue(DISPLACEMENT)

## transfer results to the model_part post
model.isogeometric_post_utility.TransferNodalResults(DISPLACEMENT, model_part_post)
model.isogeometric_post_utility.TransferNodalResults(REACTION, model_part_post)
model.isogeometric_post_utility.TransferNodalResults(STRESSES, model_part_post)

## write to gid post file
gid_io.InitializeMesh( time )
post_mesh = model_part_post.GetMesh()
print 'mesh = ', post_mesh
#gid_io.WriteNodeMesh( post_mesh )
gid_io.WriteMesh( post_mesh )
print("mesh written...")
gid_io.FinalizeMesh()
gid_io.InitializeResults( time, post_mesh )
print("write nodal displacements")
gid_io.WriteNodalResults(DISPLACEMENT, model_part_post.Nodes, time, 0)
gid_io.WriteNodalResults(REACTION, model_part_post.Nodes, time, 0)
gid_io.WriteNodalResults(STRESSES, model_part_post.Nodes, time, 0)
gid_io.FinalizeResults()


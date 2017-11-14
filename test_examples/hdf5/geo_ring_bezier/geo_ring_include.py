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
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.DiscontinuitiesApplication import *
from KratosMultiphysics.EkateAuxiliaryApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ExternalConstitutiveLawsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path ):
        #setting the domain size for the problem to be solved
        self.domain_size = 3
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("isogeometric_simulation")
        self.path = path
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
        number_of_time_steps = 1
        self.analysis_parameters = []
        # content of analysis_parameters:
        # perform_contact_analysis_flag
        # penalty value for normal contact
        # maximum number of uzawa iterations
        # friction coefficient
        # penalty value for frictional contact
        # contact_double_check_flag
        # contact_ramp_penalties_flag
        # maximum penalty value for normal contact
        # ramp criterion for normal contact
        # ramp factor for normal contact
        # maximum penalty value for frictional contact
        # ramp criterion for frictional contact
        # ramp factor for frictional contact
        # analysis type: static (0), quasi-static (1) or dynamic (2)
        perform_contact_analysis_flag = False
        penalty = 0.0
        maxuzawa = 0.0
        friction = 0.0
        frictionpenalty = 0.0
        contact_double_check_flag = False
        contact_ramp_penalties_flag = False
        maxpenalty = 0.0
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = 0.0
        fricrampcriterion = 0.0
        fricrampfactor = 0.0
        self.analysis_parameters.append(perform_contact_analysis_flag)
        self.analysis_parameters.append(penalty)
        self.analysis_parameters.append(maxuzawa)
        self.analysis_parameters.append(friction)
        self.analysis_parameters.append(frictionpenalty)
        self.analysis_parameters.append(contact_double_check_flag)
        self.analysis_parameters.append(contact_ramp_penalties_flag)
        self.analysis_parameters.append(maxpenalty)
        self.analysis_parameters.append(rampcriterion)
        self.analysis_parameters.append(rampfactor)
        self.analysis_parameters.append(fricmaxpenalty)
        self.analysis_parameters.append(fricrampcriterion)
        self.analysis_parameters.append(fricrampfactor)
        #PrintSparsityInfoFlag
        self.analysis_parameters.append(False)
        self.analysis_parameters.append(0)
        
        abs_tol =            0
        #rel_tol =            0
        rel_tol = 1e-10
        
        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        
        self.model_part.AddNodalSolutionStepVariable(LINE_LOAD)
        self.model_part.AddNodalSolutionStepVariable(STRESSES)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostAscii
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = IsogeometricModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.isogeometric_post_utility = IsogeometricClassicalPostUtility(self.model_part)
        self.generate_post_model_part = False
        self.meshWritten = True
        self.solver.CalculateReactionFlag = True
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        print "input data read OK"
        #print "+++++++++++++++++++++++++++++++++++++++"
        #for node in self.model_part.Nodes:
        #    print node
        #print "+++++++++++++++++++++++++++++++++++++++"
        
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################        
        structural_solver_advanced.AddDofs( self.model_part )
        #ekate_solver_parallel.AddDofs( self.model_part )

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = MKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2);
        
        #defined linear solver for post-processing
#        self.solver_post = SkylineLUFactorizationSolver()
        self.solver_post = MKLPardisoSolver()

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )
        
    def SetUpActivationLevels( self, model_part, activation_list, cond_activation_list ):
        for element in self.model_part.Elements:
            element.SetValue(ACTIVATION_LEVEL, activation_list[element.Id])
        for condition in self.model_part.Conditions:
            if( not (condition.GetValue(IS_TYING_MASTER) or condition.GetValue(IS_CONTACT_MASTER) ) ):
                condition.SetValue(ACTIVATION_LEVEL, activation_list[cond_activation_list[condition.Id-1]])

    def write_restart_file( self, time ):
        print("------------> restart file written for time step: "+str(time))
        self.restart_utility.ChangeFileName(problem_name+str(time))
        self.restart_utility.StoreNodalVariables(model_part)
        self.restart_utility.StoreInSituStress(model_part)
        self.restart_utility.StoreConstitutiveLawVariables(model_part)

    def WriteOutput( self, time ):
        ##transfer the results from Gauss point to nodes
        self.isogeometric_post_utility.TransferVariablesToNodes(STRESSES, self.model_part, self.solver_post)
        
        ##write result to hdf5
        hdf5_post = HDF5PostUtility(self.problem_name + "_" + str(time) + ".h5")
        hdf5_post.WriteNodes(self.model_part)
        hdf5_post.WriteNodalResults(DISPLACEMENT, self.model_part)
        hdf5_post.WriteNodalResults(REACTION, self.model_part)
        hdf5_post.WriteNodalResults(STRESSES, self.model_part)
                
    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStrainSD() )
        self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        
        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        ## ELEMENTS on layers ############################################
        
        ## NODES on layers ###############################################
        
        ## CONTACT MASTER NODES ##########################################

        ## CONTACT SLAVE NODES ###########################################

        ## INNER BOUNDARY NODES ##########################################

        ##################################################################
        print "layer sets stored"        
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        #self.SetUpActivationLevels( self.model_part, self.activation_flags, self.cond_activation_flags )
        self.deac.Initialize( self.model_part )
        self.model_part.Check(self.model_part.ProcessInfo)
        print "model successfully checked and initialized"
        
    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        if( mpi.rank == 0 ):
            self.mergefile.close()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
#        self.deac.Reactivate( self.model_part, from_reac, to_reac )
#        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
        
##################################################################

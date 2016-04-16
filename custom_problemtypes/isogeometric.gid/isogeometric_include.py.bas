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
from KratosMultiphysics import **
from KratosMultiphysics.StructuralApplication import **
from KratosMultiphysics.EkateAuxiliaryApplication import **
from KratosMultiphysics.ExternalSolversApplication import **
from KratosMultiphysics.ExternalConstitutiveLawsApplication import **
from KratosMultiphysics.IncompressibleFluidApplication import **
from KratosMultiphysics.MeshingApplication import **
from KratosMultiphysics.MKLSolversApplication import **
from KratosMultiphysics.DiscontinuitiesApplication import **
from KratosMultiphysics.IsogeometricApplication import **
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
from KratosMultiphysics.mpi import **
from KratosMultiphysics.MetisApplication import **
from KratosMultiphysics.TrilinosApplication import **
from KratosMultiphysics.MultiphaseApplication import **
*endif
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
        self.model_part_post = ModelPart("isogeometric_mesh")
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
*if(strcmp(GenData(analysis_type),"static")==0)
*set var analysistype(int)=0
*endif
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
*set var analysistype(int)=1
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
*set var analysistype(int)=2
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        perform_contact_analysis_flag = True
        # performing contact analysis: reading contact parameters
*set var penalty(real)=GenData(Penalty_Value,real)
        penalty = *penalty
*set var maxuzawa(int)=GenData(Max_Uzawa_Iterations,int)
        maxuzawa = *maxuzawa
*set var friction(real)=GenData(Friction_Coefficient,real)
        friction = *friction
*set var frictionpenalty(real)=GenData(Friction_Penalty_Value,real)
        frictionpenalty = *frictionpenalty
*if(strcmp(GenData(Bilateral_Contact),"1")==0)
        contact_double_check_flag = True
*else
        contact_double_check_flag = False
*endif
*if(strcmp(GenData(Ramp_Penalties),"1")==0)
        contact_ramp_penalties_flag = True
*set var maxpenalty(real)=GenData(Maximum_Penalty,real)
        maxpenalty = *maxpenalty
*set var rampcriterion(real)=GenData(Ramp_Criterion,real)
        rampcriterion = *rampcriterion
*set var rampfactor(real)=GenData(Ramp_Factor,real)
        rampfactor = *rampfactor
*set var fricmaxpenalty(real)=GenData(Friction_Maximum_Penalty,real)
        fricmaxpenalty = *fricmaxpenalty
*set var fricrampcriterion(real)=GenData(Friction_Ramp_Criterion,real)
        fricrampcriterion = *fricrampcriterion
*set var fricrampfactor(real)=GenData(Friction_Ramp_Factor,real)
        fricrampfactor = *fricrampfactor
*else
        contact_ramp_penalties_flag = False
        maxpenalty = penalty
        rampcriterion = 0.0
        rampfactor = 0.0
        fricmaxpenalty = penalty
        fricrampcriterion = 0.0
        fricrampfactor = 0.0
*endif
*else
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
*endif
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
*if(strcmp(GenData(Plot_Matrix_Structure),"0")==0)
        self.analysis_parameters.append(False)
*else
        self.analysis_parameters.append(True)
*endif
        self.analysis_parameters.append(*analysistype)
        
*set var abstol(real)=GenData(Absolute_Tolerance,real)
        abs_tol = *abstol
*set var reltol(real)=GenData(Relative_Tolerance,real)
        #rel_tol = *reltol
        rel_tol = 1e-10
        
        ## generating solver
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        import trilinos_contact_structural_solver
        self.solver = trilinos_contact_structural_solver.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        trilinos_contact_structural_solver.AddVariables( self.model_part )

        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = TrilinosGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        #self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        #self.gid_io = GidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        #### Reading and partitioning 
        number_of_partitions = mpi.size #we set it equal to the number of processors       
        if( mpi.rank == 0 ):
            self.mergefile = open( self.path+self.problem_name+"_merge_results.bch", 'w' )
            self.mergefile.write("Postprocess\n")
            self.mergefile.write("mescape\n \n")
        mpi.world.barrier
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        self.ReadContactNodes()
        contact_indices = IndicesVector()
        contact_indices[:] = self.contact_nodes
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        #metis_partitioning_process = MetisPartitioningProcessQuadratic(self.model_part, self.model_part_io, number_of_partitions)
        metis_partitioning_process = MetisContactPartitioningProcess(self.model_part, self.model_part_io, number_of_partitions, contact_indices)
*else
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        metis_partitioning_process = MetisPartitioningProcessQuadratic(self.model_part, self.model_part_io, number_of_partitions)
        #metis_partitioning_process = MetisPartitioningProcess(self.model_part, self.model_part_io, number_of_partitions)
*endif
        metis_partitioning_process.Execute()
        self.meshWritten = False
        mpi.world.barrier
*else
*if(strcmp(GenData(Parallel_Execution),"serial")==0)
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
*else
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
*endif

*if(strcmp(GenData(Displacements),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(DISPLACEMENT)
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(REACTION)
*endif
*if(strcmp(GenData(PK2_Stresses),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(PK2_STRESS_TENSOR)
*endif
*if(strcmp(GenData(Green_Lagrange_Strains),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(GREEN_LAGRANGE_STRAIN_TENSOR)
        self.model_part_post.AddNodalSolutionStepVariable(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
*endif
*if(strcmp(GenData(Stresses),"1")==0)
        self.model_part.AddNodalSolutionStepVariable(STRESSES)
        self.model_part_post.AddNodalSolutionStepVariable(STRESSES)
*endif
*if(strcmp(GenData(Jack_Forces),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(JACK_FORCE)
*endif
*if(strcmp(GenData(Insitu_Stress),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(PRESTRESS)
*endif
*if(strcmp(GenData(Plastic_Strains),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(PLASTIC_STRAIN_VECTOR)
        self.model_part_post.AddNodalSolutionStepVariable(PLASTICITY_INDICATOR)
*endif
*if(strcmp(GenData(Internal_Variables),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(INTERNAL_VARIABLES)
*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(AIR_PRESSURE)
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(WATER_PRESSURE)
        self.model_part_post.AddNodalSolutionStepVariable(EXCESS_PORE_WATER_PRESSURE)
*endif
*if(strcmp(GenData(Saturation),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(SATURATION)
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        self.model_part_post.AddNodalSolutionStepVariable(PARTITION_INDEX)
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        self.model_part_post.AddNodalSolutionStepVariable(CONTACT_PENETRATION)
*endif

        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
*if(strcmp(GenData(Output_Format),"ASCII")==0)
        post_mode = GiDPostMode.GiD_PostAscii
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        multi_file_flag = MultiFileFlag.MultipleFiles
*else
        multi_file_flag = MultiFileFlag.SingleFile
*endif
*else
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
*endif
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.isogeometric_post_utility = IsogeometricClassicalPostUtility(self.model_part)
        self.generate_post_model_part = False
        self.meshWritten = False
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        self.solver.CalculateReactionFlag = True
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        self.number_of_elements = len(self.model_part.Elements)
        mpi.world.barrier
        self.all_numbers_of_elements = mpi.allgather( mpi.world, self.number_of_elements )
        self.number_of_all_elements = 0
        for item in self.all_numbers_of_elements:
            self.number_of_all_elements = self.number_of_all_elements + item
        self.element_activation_levels = [0]**(self.number_of_all_elements+1)
        for element in self.model_part.Elements:
            self.element_activation_levels[element.Id] = element.GetValue(ACTIVATION_LEVEL)
        mpi.world.barrier
        for i in range(0,self.number_of_all_elements):
            levels_from_all_ranks = mpi.allgather( mpi.world, self.element_activation_levels[i] )
            for item in levels_from_all_ranks:
                if( item != 0 ):
                    self.element_activation_levels[i] = item
        mpi.world.barrier
        #print( self.element_activation_levels )
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                if( (int(val_set[1]) in self.model_part.Conditions) ):
                    self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.element_activation_levels[int(val_set[2])] )
                    print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.element_activation_levels[int(val_set[2])]) )
                    print( self.model_part.Conditions[int(val_set[1])].GetValue( ACTIVATION_LEVEL ) )
        mpi.world.barrier
        if( mpi.rank == 0 ):
            print("rank 0:")
            print(self.model_part)
        #    for condition in self.model_part.Conditions:
        #        print( condition.GetValue(ACTIVATION_LEVEL) )
        mpi.world.barrier
        if( mpi.rank == 1 ):
            print("rank 1:")
            print(self.model_part)
        #    for condition in self.model_part.Conditions:
        #        print( condition.GetValue(ACTIVATION_LEVEL) )
        mpi.world.barrier
        if( mpi.rank == 2 ):
            print("rank 2:")
            print(self.model_part)
        #    for condition in self.model_part.Conditions:
        #        print( condition.GetValue(ACTIVATION_LEVEL) )
        mpi.world.barrier
        if( mpi.rank == 3 ):
            print("rank 3:")
            print(self.model_part)
        #    for condition in self.model_part.Conditions:
        #        print( condition.GetValue(ACTIVATION_LEVEL) )

*else
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
*endif
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
*if(strcmp(GenData(Perform_MultiFlow_Analysis),"1")==0)
        for node in self.model_part.Nodes:
            node.AddDof( WATER_PRESSURE )
*if(strcmp(GenData(Perform_ThreePhase_Analysis),"1")==0)
            node.AddDof( AIR_PRESSURE )
*endif
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        trilinos_contact_structural_solver.AddDofs( self.model_part )
        #trilinos_structural_solver_static.AddDofs( self.model_part )
*else
*if(strcmp(GenData(Parallel_Execution),"serial")==0)
        structural_solver_advanced.AddDofs( self.model_part )
        #ekate_solver_parallel.AddDofs( self.model_part )
*else
        structural_solver_advanced.AddDofs( self.model_part )
*endif
*endif

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        plinear_solver = TrilinosLinearSolver()
*else
*if(strcmp(GenData(Parallel_Execution),"serial")==0)
*if(strcmp(GenData(Solver),"BiCGStabLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        preconditioner = ILU0Preconditioner()
        #preconditioner = Preconditioner()
        plinear_solver = BICGSTABSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"CGLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        preconditioner = ILU0Preconditioner()
        #preconditioner = Preconditioner()
        #plinear_solver = DeflatedCGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner,1)
        plinear_solver = CGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"SuperLU")==0)
        plinear_solver = SuperLUSolver()
*endif
*if(strcmp(GenData(Solver),"SkylineLUFactorizationSolver")==0)
        plinear_solver = SkylineLUFactorizationSolver()
*endif
*if(strcmp(GenData(Solver),"GMRESSolver")==0)
        plinear_solver = MKLGMRESSolver()
*endif
*if(strcmp(GenData(Solver),"Pardiso")==0)
        plinear_solver = MKLPardisoSolver()
*endif
*else
*if(strcmp(GenData(Solver),"SuperLU")==0)
        plinear_solver = ParallelSuperLUSolver()
*endif
*if(strcmp(GenData(Solver),"Pardiso")==0)
        plinear_solver = MKLPardisoSolver()
*endif
*if(strcmp(GenData(Solver),"BiCGStabLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        preconditioner = ILU0Preconditioner()
        #preconditioner = Preconditioner()
        plinear_solver = BICGSTABSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
*endif
*if(strcmp(GenData(Solver),"CGLinearSolver")==0)
        #preconditioner = DiagonalPreconditioner()
        preconditioner = ILU0Preconditioner()
        #preconditioner = Preconditioner()
        plinear_solver = DeflatedCGSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int),preconditioner,1)
*endif
*if(strcmp(GenData(Solver),"GMRESSolver")==0)
        #preconditioner = Preconditioner()
        #plinear_solver = GMRESSolver(*GenData(Solver_Tolerance,real), *GenData(Max_Solver_Iterations,int), preconditioner)
        plinear_solver = ParallelMKLGMRESSolver()
*endif
*if(strcmp(GenData(Solver),"SkylineLUFactorizationSolver")==0)
        plinear_solver = SkylineLUFactorizationSolver()
*endif
*endif
        self.solver.structure_linear_solver = plinear_solver
*endif
        self.solver_post = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2);

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

#    def restart_time_step( self, time, Dt ):
#        print("############ time step solution has to be restarted ############")
#        time = time-Dt
#        model_part.CloneTimeStep(time)
#        for step in range(1,11):
#            time = time+ Dt/10.0
#            model_part.CloneTimeStep(time)
#            #####################################################################################################
#            model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
#            model_part.ProcessInfo.SetValue( FIRST_TIME_STEP, False )
#            #####################################################################################################
#            solver.Solve()
#            print("~~~~~~~~~~~~~~ RESTARTED STEP ( DT= "+str(Dt/10.0)+" / Step= "+str(step)+" ) ~~~~~~~~~~~~~~")
#        print("############ restart finished ############")
#
#    def write_to_file( self, time ):
#        for i in range(0, len(self.layer_nodes_sets['top'])):
#            settlements.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetSolutionStepValue(DISPLACEMENT_Z))+"\n")
#    for i in range(0, len(layer_nodes_sets['side'])):
#        pressure_air.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(AIR_PRESSURE))+"\n")
#        pressure_water.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(WATER_PRESSURE))+"\n")
#
            
    def WriteMaterialParameters( self, time, indices ):
        self.gid_io.OpenResultFile( self.path+self.problem_name, GiDPostMode.GiD_PostBinary)
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        #self.gid_io.ChangeOutputName( self.path+self.problem_name +str(time), GiDPostMode.GiD_PostBinary )
*endif
        for index in indices:
            self.gid_io.SuperPrintOnGaussPoints(MATERIAL_PARAMETERS, self.model_part_post, time, index)
        self.gid_io.CloseResultFile()

    def WriteMonitoringSectionResults( self, time ):
        outfile = open("step_"+str(time)+".dat",'w')
        outfile.write("ekate result file for step "+str(time)+"\n")
        outfile.close()
        
    def WriteOutput( self, time ):
        if self.generate_post_model_part == False:
*if(strcmp(GenData(Post_Element_Type),"Triangle")==0)
            post_element_type = PostElementType.Triangle
*endif
*if(strcmp(GenData(Post_Element_Type),"Quadrilateral")==0)
            post_element_type = PostElementType.Quadrilateral
*endif
*if(strcmp(GenData(Post_Element_Type),"Tetrahedra")==0)
            post_element_type = PostElementType.Tetrahedra
*endif
*if(strcmp(GenData(Post_Element_Type),"Hexahedra")==0)
            post_element_type = PostElementType.Hexahedra
*endif
            self.isogeometric_post_utility.GenerateModelPart(self.model_part_post, post_element_type)
            self.generate_post_model_part = True
            print 'Generate PostModelPart completed'
        if self.generate_post_model_part == True:
*if(strcmp(GenData(Displacements),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(DISPLACEMENT, self.model_part_post)
*endif
*if(strcmp(GenData(Reactions),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(REACTION, self.model_part_post)
*endif
*if(strcmp(GenData(PK2_Stresses),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(PK2_STRESS_TENSOR, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Green_Lagrange_Strains),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part_post, self.solver_post)
            self.isogeometric_post_utility.TransferIntegrationPointResults(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Stresses),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(STRESSES, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Jack_Forces),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(JACK_FORCE, self.model_part_post)
*endif
*if(strcmp(GenData(Insitu_Stress),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(PRESTRESS, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Plastic_Strains),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(PLASTIC_STRAIN_VECTOR, self.model_part_post, self.solver_post)
            self.isogeometric_post_utility.TransferIntegrationPointResults(PLASTICITY_INDICATOR, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Internal_Variables),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(INTERNAL_VARIABLES, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(AIR_PRESSURE, self.model_part_post)
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(WATER_PRESSURE, self.model_part_post)
            self.isogeometric_post_utility.TransferIntegrationPointResults(EXCESS_PORE_WATER_PRESSURE, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Saturation),"1")==0)
            self.isogeometric_post_utility.TransferIntegrationPointResults(SATURATION, self.model_part_post, self.solver_post)
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
            self.isogeometric_post_utility.TransferNodalResults(PARTITION_INDEX, self.model_part_post)
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
            self.isogeometric_post_utility.TransferNodalResults(CONTACT_PENETRATION, self.model_part_post)
*endif
            print 'Transfer Nodal Results completed'
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        meshname = time+0.001**(mpi.rank+1)
        self.gid_io.InitializeMesh( meshname )
*else
        self.gid_io.InitializeMesh( time )
*endif
        mesh = self.model_part_post.GetMesh()
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
*else
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        if( self.meshWritten == False ):
            meshname = time+0.001**mpi.rank
            self.gid_io.InitializeMesh( meshname )
*else
        if( self.meshWritten == False ):
            self.gid_io.InitializeMesh( 0.0 )
*endif
            mesh = self.model_part_post.GetMesh()
            self.gid_io.WriteMesh( mesh )
            self.meshWritten = True
            self.gid_io.FinalizeMesh()
*endif
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        self.gid_io.InitializeResults( time, self.model_part_post.GetMesh() )
*else
        self.gid_io.InitializeResults( time, self.model_part_post.GetMesh() )
*endif
*else
        self.gid_io.InitializeResults( 0.0, self.model_part_post.GetMesh() )
*endif
*if(strcmp(GenData(Displacements),"1")==0)
        print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(PK2_Stresses),"1")==0)
        self.gid_io.WriteNodalResults(PK2_STRESS_TENSOR, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Green_Lagrange_Strains),"1")==0)
        self.gid_io.WriteNodalResults(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Stresses),"1")==0)
        self.gid_io.WriteNodalResults(STRESSES, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Jack_Forces),"1")==0)
        self.gid_io.WriteNodalResults(JACK_FORCE, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Insitu_Stress),"1")==0)
        self.gid_io.WriteNodalResults(PRESTRESS, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Plastic_Strains),"1")==0)
        self.gid_io.WriteNodalResults(PLASTIC_STRAIN_VECTOR, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(PLASTICITY_INDICATOR, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Internal_Variables),"1")==0)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 1)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 2)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 3)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 4)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 5)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 6)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 7)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 8)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 9)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 10)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 11)
        self.gid_io.WriteNodalResults(INTERNAL_VARIABLES, self.model_part_post.Nodes, time, 12)

*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
        self.gid_io.WriteNodalResults(AIR_PRESSURE, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
        self.gid_io.WriteNodalResults(WATER_PRESSURE, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(EXCESS_PORE_WATER_PRESSURE, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Saturation),"1")==0)
        self.gid_io.WriteNodalResults(SATURATION, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        self.gid_io.WriteNodalResults(PARTITION_INDEX, self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(Perform_Contact_Analysis),"1")==0)
        self.gid_io.WriteNodalResults(CONTACT_PENETRATION,self.model_part_post.Nodes, time, 0)
*endif
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.FinalizeResults()
*endif
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        if( mpi.rank == 0 ):
            meshname = time+0.001
            self.mergefile.write("Files Read "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
            self.mergefile.write("mescape\n")
            for rank in range(2, mpi.size+1 ):
                meshname = time+0.001**(rank)
                self.mergefile.write("Files Add "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
                self.mergefile.write("mescape\n")
            self.mergefile.write("Files SaveAll BinMeshesSets\n")
            self.mergefile.write(self.path+"merged_"+self.problem_name+"_"+ str(time)+".bin\n")
            self.mergefile.write("mescape\n")
            self.mergefile.write("\n")
*endif
                
    def InitializeModel( self, model_data ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        
        #TODO: handling the properties
        
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
*if(strcmp(GenData(Parallel_Execution),"mpi")==0)
        self.deac = TrilinosDeactivationUtility()
*else
        self.deac = DeactivationUtility()
*endif
        #self.SetUpActivationLevels( self.model_part, self.activation_flags, self.cond_activation_flags )
        self.deac.Initialize( self.model_part ) #elements and conditions is initialized here
        self.model_part.Check(self.model_part.ProcessInfo) #important to check the negativeness of Jacobian
        print "model successfully checked and initialized"
        
    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        if( mpi.rank == 0 ):
            self.mergefile.close()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
        
##################################################################

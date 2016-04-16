//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 18, 2013 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/legacy_structural_app_vars.h"
#include "includes/element.h"
#include "includes/condition.h"

#include "custom_elements/kinematic_linear_nurbs.h"
#include "custom_elements/kinematic_linear_isogeometric.h"
#include "custom_elements/total_lagrangian_isogeometric.h"
#include "custom_elements/unsaturated_soils_element_2phase_small_strain_isogeometric.h"
#include "custom_conditions/line_force_isogeometric.h"
#include "custom_conditions/line_force_isogeometric_2d.h"
#include "custom_conditions/face_load_isogeometric.h"
#include "custom_conditions/face_pressure_isogeometric.h"
#include "custom_conditions/slave_contact_face_3D_isogeometric.h"
#include "custom_conditions/master_contact_face_3D_isogeometric.h"
#include "custom_geometries/geo_1d_nurbs.h"
#include "custom_geometries/geo_2d_nurbs.h"
#include "custom_geometries/geo_3d_nurbs.h"
#include "custom_geometries/geo_1d_bezier.h"
#include "custom_geometries/geo_2d_bezier.h"
#include "custom_geometries/geo_2d_bezier_3.h"
#include "custom_geometries/geo_3d_bezier.h"
#include "structural_application/custom_elements/kinematic_linear.h"
#include "structural_application/custom_elements/unsaturated_soils_element_2phase_small_strain.h"
#include "phase_field_application/custom_elements/phase_field_fracture.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_WEIGHT )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_1 )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_2 )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_3 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_1 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_2 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_3 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_1 ) //number of control points in 1st direction
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_2 ) //number of control points in 2nd direction
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_3 ) //number of control points in 3rd direction
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_1 ) //number of mesh points along 1st direction in post-processing
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_2 ) //number of mesh points along 2nd direction in post-processing
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_3 ) //number of mesh points along 3rd direction in post-processing
    KRATOS_DEFINE_VARIABLE( Matrix, EXTRACTION_OPERATOR )
    KRATOS_DEFINE_VARIABLE( Matrix, EXTRACTION_OPERATOR_MCSR )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_ROWPTR )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_COLIND )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_VALUES )
    
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_COORDINATES)

    KRATOS_DEFINE_VARIABLE( double, STABILISATION_FACTOR )
    KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS )
    
	///@} 
	///@name Type Definitions
	///@{ 

	///@} 
	///@name  Enum's
	///@{

	///@}
	///@name  Functions 
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosIsogeometricApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosDiscontinuitiesApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosIsogeometricApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosIsogeometricApplication();

		/// Destructor.
		virtual ~KratosIsogeometricApplication()
		{}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



		///@}
		///@name Access
		///@{ 


		///@}
		///@name Inquiry
		///@{


		///@}      
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosIsogeometricApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in KratosIsogeometricApplication");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
      }


		///@}      
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


		///@} 
		///@name Protected  Access 
		///@{ 


		///@}      
		///@name Protected Inquiry 
		///@{ 


		///@}    
		///@name Protected LifeCycle 
		///@{ 


		///@}

	private:
		///@name Static Member Variables 
		///@{ 


		///@} 
		///@name Member Variables 
		///@{ 
        
//        const KinematicLinearNURBS mKinematicLinearGeo1dNURBS; //TODO: doesn't work, segmentation fault error
//        const KinematicLinearNURBS mKinematicLinearGeo2dNURBS;
//        const KinematicLinearNURBS mKinematicLinearGeo3dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo1dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo2dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo3dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo1dBezier;
        const KinematicLinearIsogeometric mKinematicLinearGeo2dBezier;
        const KinematicLinearIsogeometric mKinematicLinearGeo3dBezier;
        const TotalLagrangianIsogeometric mTotalLagrangianGeo3dBezier;
        const UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier;
        
        const LineForceIsogeometric mLineLoadNURBS;
        const LineForceIsogeometric2D mLineLoadNURBS2D;
        const LineForceIsogeometric mLineLoadBezier;
        const LineForceIsogeometric2D mLineLoadBezier2D;

        const FaceLoadIsogeometric mFaceLoadNURBS;
        const FaceLoadIsogeometric mFaceLoadBezier;
        const FacePressureIsogeometric mFacePressureNURBS;
        const FacePressureIsogeometric mFacePressureBezier;
        
        const MasterContactFace3DIsogeometric mMasterContactFace3DBezier;
        const SlaveContactFace3DIsogeometric mSlaveContactFace3DBezier;

        const KinematicLinear mKinematicLinearBezier2D;
        const KinematicLinear mKinematicLinearBezier3D;
        const UnsaturatedSoilsElement_2phase_SmallStrain mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D;

        const PhaseFieldFracture mPhaseFieldFractureBezier2D;
        const PhaseFieldFracture mPhaseFieldFractureBezier3D;

		///@} 
		///@name Private Operators
		///@{ 


		///@} 
		///@name Private Operations
		///@{ 


		///@} 
		///@name Private  Access 
		///@{ 


		///@}    
		///@name Private Inquiry 
		///@{ 


		///@}    
		///@name Un accessible methods 
		///@{ 

		/// Assignment operator.
		KratosIsogeometricApplication& operator=(KratosIsogeometricApplication const& rOther);

		/// Copy constructor.
		KratosIsogeometricApplication(KratosIsogeometricApplication const& rOther);


		///@}    

	}; // Class KratosIsogeometricApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED  defined 



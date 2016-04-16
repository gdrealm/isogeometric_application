//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2013 $
//   Revision:            $Revision: 1.1 $
//
// 


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "isogeometric_application.h"
#include "custom_geometries/geo_1d_bezier.h"
#include "custom_geometries/geo_2d_bezier.h"
#include "custom_geometries/geo_2d_bezier_3.h"
#include "custom_geometries/geo_3d_bezier.h"


namespace Kratos
{
    
    KRATOS_CREATE_VARIABLE( Vector, NURBS_WEIGHT )
    KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_1 )
    KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_2 )
    KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_3 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_1 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_2 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_3 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_1 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_2 )
    KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_3 )
    KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_1 )
    KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_2 )
    KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_3 )
    KRATOS_CREATE_VARIABLE( Matrix, EXTRACTION_OPERATOR )
    KRATOS_CREATE_VARIABLE( Matrix, EXTRACTION_OPERATOR_MCSR )
    KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_ROWPTR )
    KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_COLIND )
    KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_VALUES )

    KRATOS_CREATE_VARIABLE( double, STABILISATION_FACTOR )
    KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS )

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )

    KratosIsogeometricApplication::KratosIsogeometricApplication()
    :
      mKinematicLinearGeo1dNURBS( 0, Element::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo2dNURBS( 0, Element::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo3dNURBS( 0, Element::GeometryType::Pointer( new Geo3dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo1dBezier( 0, Element::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mKinematicLinearGeo2dBezier( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mKinematicLinearGeo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mTotalLagrangianGeo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mLineLoadNURBS( 0, Condition::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mLineLoadNURBS2D( 0, Condition::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mLineLoadBezier( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mLineLoadBezier2D( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mFaceLoadNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFaceLoadBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mFacePressureNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFacePressureBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mMasterContactFace3DBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mSlaveContactFace3DBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKinematicLinearBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mKinematicLinearBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mPhaseFieldFractureBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mPhaseFieldFractureBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    {}
    
 	void KratosIsogeometricApplication::Register()
 	{
 		// calling base class to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosIsogeometricApplication... " << std::endl;
        
        // register elements
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo1dNURBS", mKinematicLinearGeo1dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo2dNURBS", mKinematicLinearGeo2dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo3dNURBS", mKinematicLinearGeo3dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo1dBezier", mKinematicLinearGeo1dBezier )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo2dBezier", mKinematicLinearGeo2dBezier )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo3dBezier", mKinematicLinearGeo3dBezier )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianGeo3dBezier", mTotalLagrangianGeo3dBezier )
        KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier", mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier )
        
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS3D", mLineLoadNURBS )
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS2D", mLineLoadNURBS2D )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier3D", mLineLoadBezier )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier2D", mLineLoadBezier2D )
        KRATOS_REGISTER_CONDITION( "FaceLoadNURBS", mFaceLoadNURBS )
        KRATOS_REGISTER_CONDITION( "FaceLoadBezier", mFaceLoadBezier )
        KRATOS_REGISTER_CONDITION( "FacePressureNURBS", mFacePressureNURBS )
        KRATOS_REGISTER_CONDITION( "FacePressureBezier", mFacePressureBezier )
        KRATOS_REGISTER_CONDITION( "MasterContactFace3DBezier", mMasterContactFace3DBezier )
        KRATOS_REGISTER_CONDITION( "SlaveContactFace3DBezier", mSlaveContactFace3DBezier )

        KRATOS_REGISTER_ELEMENT( "KinematicLinearBezier2D", mKinematicLinearBezier2D )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearBezier3D", mKinematicLinearBezier3D )
        KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement_2phase_SmallStrainBezier3D", mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D )
        KRATOS_REGISTER_ELEMENT( "PhaseFieldFractureBezier2D", mPhaseFieldFractureBezier2D )
        KRATOS_REGISTER_ELEMENT( "PhaseFieldFractureBezier3D", mPhaseFieldFractureBezier3D )

        //register variables
        KRATOS_REGISTER_VARIABLE( NURBS_WEIGHT )
        KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_1 )
        KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_2 )
        KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_3 )
        KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_1 )
        KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_2 )
        KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_3 )
        KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_1 )
        KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_2 )
        KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_3 )
        KRATOS_REGISTER_VARIABLE( NUM_DIVISION_1 )
        KRATOS_REGISTER_VARIABLE( NUM_DIVISION_2 )
        KRATOS_REGISTER_VARIABLE( NUM_DIVISION_3 )
        KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR )
        KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_MCSR )
        KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_ROWPTR )
        KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_COLIND )
        KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_VALUES )

        KRATOS_REGISTER_VARIABLE( STABILISATION_FACTOR)
        KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS)

        // to make sure the variable imported from other application is registerered
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )

        // register the geometries
        Geo1dBezier<Node<3> > Geo1dBezierPrototype;
        Serializer::Register( "Geo1dBezier", Geo1dBezierPrototype );

        Geo2dBezier<Node<3> > Geo2dBezierPrototype;
        Serializer::Register( "Geo2dBezier", Geo2dBezierPrototype );

        Geo2dBezier3<Node<3> > Geo2dBezier3Prototype;
        Serializer::Register( "Geo2dBezier3", Geo2dBezier3Prototype );

        Geo3dBezier<Node<3> > Geo3dBezierPrototype;
        Serializer::Register( "Geo3dBezier", Geo3dBezierPrototype );
 	}

}  // namespace Kratos.



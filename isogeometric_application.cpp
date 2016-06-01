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

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )

    KratosIsogeometricApplication::KratosIsogeometricApplication()
    {}

 	void KratosIsogeometricApplication::Register()
 	{
 		// calling base class to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosIsogeometricApplication... " << std::endl;
        
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


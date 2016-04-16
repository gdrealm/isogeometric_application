#pragma once

//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Simple_Homogeneous.h>
//#include <CGAL/MP_Float.h>
//#include <CGAL/Quotient.h>
//#include <CGAL/basic.h>
//#include <CGAL/Arr_segment_traits_2.h>
//#include <CGAL/Arrangement_2.h>
//#include <CGAL/Arr_extended_dcel.h>
//#include <CGAL/Arr_naive_point_location.h>
//#include <CGAL/Arr_dcel_base.h>
//#include <CGAL/squared_distance_3.h>
//#include <CGAL/Arr_naive_point_location.h>
//#include <CGAL/Arr_trapezoid_ric_point_location.h>
//#include <CGAL/Arr_landmarks_point_location.h>
//#include <CGAL/Arr_accessor.h>

namespace geomodcore {

    /** 
    Provides exact geometric predicates, but geometric constructions may be inexact due to round-off errors. 
    It is however enough for most Cgal algorithms, and faster than both Exact_predicates_exact_constructions_kernel
    and Exact_predicates_exact_constructions_kernel_with_sqrt. 
    */
    //typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
    //typedef Kernel::FT kdouble;        // Kernel-double-Datentyp

    //typedef CGAL::Quotient<CGAL::MP_Float> FieldNumberType;        // Numerischer Fieldnumber-Typ
    //typedef CGAL::MP_Float gdouble;                            // Numerischer Ringnumber-Typ
    //    typedef CGAL::Simple_homogeneous<double> HKernel;    // Kernel-Typ 

    const unsigned int GDimension = 3;    
    const unsigned int PDimension = 2;
    const unsigned int TSplineDegree = 3;            // kubische T-Splines global vorgeben, dies ist Template-Parameter

    typedef boost::array<double, GDimension>                        GPoint;                            // nD-Geometrie-Punkte
    typedef boost::array<double, PDimension>                        PPoint;                            // nD-Parameter/Bild-Punkte
    typedef boost::array<double, PDimension>                        IPoint;                            // nD-Parameter/Bild-Punkte
    typedef boost::array<double, 2>                                 PPoint2d;                        // 2D-Parameter/Bild-Punkte
    typedef boost::array<double, 3>                                 PPoint3d;                        // 3D-Parameter/Bild-Punkte

    typedef boost::array<double, ((GDimension+1)*GDimension)>       AffTransformG;                    // nxn Matrix
    typedef boost::array<double, ((PDimension+1)*PDimension)>       AffTransformP;                    // nxn Matrix
    typedef boost::array<double, ((PDimension+1)*PDimension)>       AffTransformI;                    // nxn Matrix

    typedef boost::array<double, 7>                                 AffTransform2d;                    // 3x3 Matrix mit m20 - m21 = 0, m22 = 1 für Parameter-Transformation
    typedef boost::array<double, 13>                                AffTransform3d;                    // 4x4 Matrix mit m30 - m32 = 0, m33 = 1 für Parameter-Transformation

    typedef std::pair<double,double>                                MinMaxP1d;
    
    typedef std::pair<IPoint,IPoint>                                MinMaxI;                        // Paar mininmaler/maximaler Bildkoordinaten
    typedef std::pair<PPoint,PPoint>                                MinMaxP;                        // Paar mininmaler/maximaler Bildkoordinaten
    typedef std::pair<GPoint,GPoint>                                MinMaxG;                        // Paar mininmaler/maximaler Geometriekoordinaten

    typedef std::vector<double>                                     Knots;
    typedef boost::array<Knots, PDimension>                         KnotsN;                            // n-dimensionaler Knotenvektor
    
    typedef boost::array<double,TSplineDegree+2>                    TKnots;
    typedef boost::array<TKnots, PDimension>                        TKnotsN;                        // n-dimensionaler Knotenvektor für T-Splines

    GPoint transform(const GPoint& p, const AffTransformG& T)
    {
        GPoint tmp = p;
        for (int j = 0; j < GDimension - 1; j++)
        {
            for (int i = 0; i < GDimension; i++)
            {
                tmp[j] += p[i] * T[j * GDimension + i];
            }
        }
        return tmp;
    };

    PPoint transform(const PPoint& p, const AffTransformP& T)
    {
        PPoint tmp = p;
        for (int j = 0; j < PDimension - 1; j++)
        {
            for (int i = 0; i < PDimension; i++)
            {
                tmp[j] += p[i] * T[j * PDimension + i];
            }
        }
        return tmp;
    };
}

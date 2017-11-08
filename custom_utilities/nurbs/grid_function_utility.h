//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/patch.h"
#include "custom_utilities/nurbs/nurbs_patch.h"

namespace Kratos
{

/**
This class is a library to generate various grid function for typical computational mechanics problems.
 */
class GridFunctionUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunctionUtility);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;

    /// Default constructor
    GridFunctionUtility() {}

    /// Destructor
    virtual ~GridFunctionUtility() {}

    /// Generate the regular equidistant control point grid. All the point has unit weight.
    template<std::size_t TDim>
    static typename GridFunction<TDim, ControlPointType>::Pointer CreateRegularControlPointGrid(
            const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& spacing)
    {
        if (TDim == 1)
        {
            typename GridFunction<1, ControlPointType>::Pointer pGrid
                = typename GridFunction<1, ControlPointType>::Pointer(new GridFunction<1, ControlPointType>(ngrid[0]+1));

            for (std::size_t i = 0; i < ngrid[0]+1; ++i)
            {
                ControlPointType Point;
                Point.SetCoordinates(start[0] + i*spacing[0], 0.0, 0.0, 1.0);
                pGrid->SetValue(i, Point);
            }
        }
        else if (TDim == 2)
        {
            typename GridFunction<2, ControlPointType>::Pointer pGrid
                = typename GridFunction<2, ControlPointType>::Pointer(new GridFunction<2, ControlPointType>(ngrid[0]+1, ngrid[1]+1));

            for (std::size_t i = 0; i < ngrid[0]+1; ++i)
            {
                for (std::size_t j = 0; j < ngrid[1]+1; ++j)
                {
                    ControlPointType Point;
                    Point.SetCoordinates(start[0] + i*spacing[0], start[1] + j*spacing[1], 0.0, 1.0);
                    pGrid->SetValue(i, j, Point);
                }
            }
        }
        else if (TDim == 3)
        {
            typename GridFunction<3, ControlPointType>::Pointer pGrid
                = typename GridFunction<3, ControlPointType>::Pointer(new GridFunction<3, ControlPointType>(ngrid[0]+1, ngrid[1]+1, ngrid[2]+1));

            for (std::size_t i = 0; i < ngrid[0]+1; ++i)
            {
                for (std::size_t j = 0; j < ngrid[1]+1; ++j)
                {
                    for (std::size_t k = 0; k < ngrid[2]+1; ++k)
                    {
                        ControlPointType Point;
                        Point.SetCoordinates(start[0] + i*spacing[0], start[1] + j*spacing[1], start[2] + k*spacing[2], 1.0);
                        pGrid->SetValue(i, j, k, Point);
                    }
                }
            }
        }

        return pGrid;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GridFunctionUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const GridFunctionUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_UTILITY_H_INCLUDED defined


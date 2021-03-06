//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/nurbs/structured_control_grid.h"


namespace Kratos
{

template<int TDim, typename TDataType>
struct ControlGridLibrary_Helper
{
    /// Generate regular control grid with a specific data type
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid);
};

/**
Helper library to generate control grid for isogeometric analysis
 */
class ControlGridLibrary
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlGridLibrary);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;

    /// Default constructor
    ControlGridLibrary() {}

    /// Destructor
    virtual ~ControlGridLibrary() {}



    /// Generate the regular equidistant control point grid. All the point has unit weight.
    template<int TDim>
    static ControlGrid<ControlPointType>::Pointer CreateStructuredControlPointGrid(
            const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& end)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }



    /// Generate the regular control point grid based on starting point and the director vector in each direction. All the point has unit weight.
    template<int TDim>
    static ControlGrid<ControlPointType>::Pointer CreateStructuredControlPointGrid(
            const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<std::vector<double> >& spacing_vectors)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }



    /// Generate regular control grid with a specific data type
    template<int TDim, typename TDataType>
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid)
    {
        return ControlGridLibrary_Helper<TDim, TDataType>::CreateStructuredZeroControlGrid(Name, ngrid);
    }



    /// Generate regular control grid with variable
    template<int TDim, class TVariableType>
    static typename ControlGrid<typename TVariableType::Type>::Pointer CreateStructuredZeroControlGrid(const TVariableType& rVariable, const std::vector<std::size_t>& ngrid)
    {
        return CreateStructuredZeroControlGrid<TDim, typename TVariableType::Type>(rVariable.Name(), ngrid);
    }



    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlGridLibrary";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<1>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& end)
{
    StructuredControlGrid<1, ControlPointType>::Pointer pGrid = StructuredControlGrid<1, ControlPointType>::Create(ngrid);

    pGrid->SetName("CONTROL_POINT");

    std::vector<double> spacing = {(end[0] - start[0]) / (ngrid[0]-1),
                                (end[1] - start[1]) / (ngrid[0]-1),
                                (end[2] - start[2]) / (ngrid[0]-1)};

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        ControlPointType Point;
        Point.SetCoordinates(start[0] + i*spacing[0], start[1] + i*spacing[1], start[2] + i*spacing[2], 1.0);
        pGrid->SetValue(i, Point);
    }

    return pGrid;
}

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<2>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& end)
{
    StructuredControlGrid<2, ControlPointType>::Pointer pGrid = StructuredControlGrid<2, ControlPointType>::Create(ngrid);

    pGrid->SetName("CONTROL_POINT");

    std::vector<double> spacing = {(end[0] - start[0]) / (ngrid[0]-1), (end[1] - start[1]) / (ngrid[1]-1)};

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        for (std::size_t j = 0; j < ngrid[1]; ++j)
        {
            ControlPointType Point;
            Point.SetCoordinates(start[0] + i*spacing[0], start[1] + j*spacing[1], 0.0, 1.0);
            pGrid->SetValue(i, j, Point);
        }
    }

    return pGrid;
}

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<3>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& end)
{
    StructuredControlGrid<3, ControlPointType>::Pointer pGrid = StructuredControlGrid<3, ControlPointType>::Create(ngrid);

    pGrid->SetName("CONTROL_POINT");

    std::vector<double> spacing = {(end[0] - start[0]) / (ngrid[0]-1),
                                (end[1] - start[1]) / (ngrid[1]-1),
                                (end[2] - start[2]) / (ngrid[2]-1)};

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        for (std::size_t j = 0; j < ngrid[1]; ++j)
        {
            for (std::size_t k = 0; k < ngrid[2]; ++k)
            {
                ControlPointType Point;
                Point.SetCoordinates(start[0] + i*spacing[0], start[1] + j*spacing[1], start[2] + k*spacing[2], 1.0);
                pGrid->SetValue(i, j, k, Point);
            }
        }
    }

    return pGrid;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<1>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<std::vector<double> >& spacing_vectors)
{
    // size check
    assert(start.size() == 3);
    assert(ngrid.size() == 1);
    assert(spacing_vectors.size() > 0);
    assert(spacing_vectors[0].size() == 3);

    typename StructuredControlGrid<1, ControlPointType>::Pointer pGrid
        = typename StructuredControlGrid<1, ControlPointType>::Pointer(new StructuredControlGrid<1, ControlPointType>(ngrid[0]));

    pGrid->SetName("CONTROL_POINT");

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        ControlPointType Point;
        Point.SetCoordinates(start[0] + i*spacing_vectors[0][0], start[1] + i*spacing_vectors[0][1], start[2] + i*spacing_vectors[0][2], 1.0);
        pGrid->SetValue(i, Point);
    }

    return pGrid;
}

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<2>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<std::vector<double> >& spacing_vectors)
{
    // size check
    assert(start.size() == 3);
    assert(ngrid.size() == 2);
    assert(spacing_vectors.size() > 1);
    assert(spacing_vectors[0].size() == 3);
    assert(spacing_vectors[1].size() == 3);

    typename StructuredControlGrid<2, ControlPointType>::Pointer pGrid
        = typename StructuredControlGrid<2, ControlPointType>::Pointer(new StructuredControlGrid<2, ControlPointType>(ngrid[0], ngrid[1]));

    pGrid->SetName("CONTROL_POINT");

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        for (std::size_t j = 0; j < ngrid[1]; ++j)
        {
            ControlPointType Point;
            double x = start[0] + i*spacing_vectors[0][0] + j*spacing_vectors[1][0];
            double y = start[1] + i*spacing_vectors[0][1] + j*spacing_vectors[1][1];
            double z = start[2] + i*spacing_vectors[0][2] + j*spacing_vectors[1][2];
            Point.SetCoordinates(x, y, z, 1.0);
            pGrid->SetValue(i, j, Point);
        }
    }

    return pGrid;
}

template<>
ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary::CreateStructuredControlPointGrid<3>(
        const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<std::vector<double> >& spacing_vectors)
{
    // size check
    assert(start.size() == 3);
    assert(ngrid.size() == 3);
    assert(spacing_vectors.size() > 2);
    assert(spacing_vectors[0].size() == 3);
    assert(spacing_vectors[1].size() == 3);
    assert(spacing_vectors[2].size() == 3);

    typename StructuredControlGrid<3, ControlPointType>::Pointer pGrid
        = typename StructuredControlGrid<3, ControlPointType>::Pointer(new StructuredControlGrid<3, ControlPointType>(ngrid[0], ngrid[1], ngrid[2]));

    pGrid->SetName("CONTROL_POINT");

    for (std::size_t i = 0; i < ngrid[0]; ++i)
    {
        for (std::size_t j = 0; j < ngrid[1]; ++j)
        {
            for (std::size_t k = 0; k < ngrid[2]; ++k)
            {
                ControlPointType Point;
                double x = start[0] + i*spacing_vectors[0][0] + j*spacing_vectors[1][0] + k*spacing_vectors[2][0];
                double y = start[1] + i*spacing_vectors[0][1] + j*spacing_vectors[1][1] + k*spacing_vectors[2][1];
                double z = start[2] + i*spacing_vectors[0][2] + j*spacing_vectors[1][2] + k*spacing_vectors[2][2];
                Point.SetCoordinates(x, y, z, 1.0);
                pGrid->SetValue(i, j, k, Point);
            }
        }
    }

    return pGrid;
}

////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TDataType>
struct ControlGridLibrary_Helper<1, TDataType>
{
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid)
    {
        typename StructuredControlGrid<1, TDataType>::Pointer pGrid
            = typename StructuredControlGrid<1, TDataType>::Pointer(new StructuredControlGrid<1, TDataType>(ngrid[0]));

        pGrid->SetName(Name);

        for (std::size_t i = 0; i < ngrid[0]; ++i)
        {
            TDataType Point(0.0);
            pGrid->SetValue(i, Point);
        }

        return pGrid;
    }
};

template<typename TDataType>
struct ControlGridLibrary_Helper<2, TDataType>
{
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid)
    {
        typename StructuredControlGrid<2, TDataType>::Pointer pGrid
            = typename StructuredControlGrid<2, TDataType>::Pointer(new StructuredControlGrid<2, TDataType>(ngrid[0], ngrid[1]));

        pGrid->SetName(Name);

        for (std::size_t i = 0; i < ngrid[0]; ++i)
        {
            for (std::size_t j = 0; j < ngrid[1]; ++j)
            {
                TDataType Point(0.0);
                pGrid->SetValue(i, j, Point);
            }
        }

        return pGrid;
    }
};

template<typename TDataType>
struct ControlGridLibrary_Helper<3, TDataType>
{
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid)
    {
        typename StructuredControlGrid<3, TDataType>::Pointer pGrid
            = typename StructuredControlGrid<3, TDataType>::Pointer(new StructuredControlGrid<3, TDataType>(ngrid[0], ngrid[1], ngrid[2]));

        pGrid->SetName(Name);

        for (std::size_t i = 0; i < ngrid[0]; ++i)
        {
            for (std::size_t j = 0; j < ngrid[1]; ++j)
            {
                for (std::size_t k = 0; k < ngrid[2]; ++k)
                {
                    TDataType Point(0.0);
                    pGrid->SetValue(i, j, k, Point);
                }
            }
        }

        return pGrid;
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const ControlGridLibrary& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED defined


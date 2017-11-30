//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include "boost/algorithm/string.hpp"

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/patch.h"

namespace Kratos
{

/**
THis class supports some operations on B-Splines patch
 */
class BSplinesPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchUtility);

    /// Type definition

    enum ReadMode
    {
        _NO_READ          = 0,
        _READ_PATCH       = 1,
        _READ_ORDER       = 2,
        _READ_NUMBER      = 3,
        _READ_KNOTS       = 4,
        _READ_COORDINATES = 5,
        _READ_WEIGHTS     = 6
    };

    /// Default constructor
    BSplinesPatchUtility() {}

    /// Destructor
    virtual ~BSplinesPatchUtility() {}

    /// Construct a higher dimension patch by connecting two patches with the straight B-Splines curve. The order of the connection curve is 1.
    /// To have higher order of the connection one needs to elevate the degree.
    /// Right now, the two sub-patches must have same parameters (knot vectors) and are B-Splines.
    template<int TDim>
    static typename Patch<TDim>::Pointer CreateConnectedPatch(typename Patch<TDim-1>::Pointer pPatch1, typename Patch<TDim-1>::Pointer pPatch2)
    {
        // check prerequisites
        if (pPatch1->pFESpace()->Type() != BSplinesFESpace<TDim-1>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, "Patch 1 is not B-Splines patch", "")

        if (pPatch2->pFESpace()->Type() != BSplinesFESpace<TDim-1>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, "Patch 2 is not B-Splines patch", "")

        // check the FESpace
        if ( !pPatch1->pFESpace()->IsCompatible(*(pPatch2->pFESpace())) )
        {
            KRATOS_THROW_ERROR(std::logic_error, "The two patches are not compatible", "")
        }

        // create the new FESpace
        typename BSplinesFESpace<TDim-1>::Pointer pFESpace1 = boost::dynamic_pointer_cast<BSplinesFESpace<TDim-1> >(pPatch1->pFESpace());
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
        for (std::size_t dim = 0; dim < TDim-1; ++dim)
        {
            pNewFESpace->SetKnotVector(dim, pFESpace1->KnotVector(dim));
            pNewFESpace->SetInfo(dim, pFESpace1->Number(dim), pFESpace1->Order(dim));
        }

        std::size_t connect_order = 1;
        typename BSplinesFESpace<TDim>::knot_container_t new_knot_vector = BSplinesFESpaceLibrary::CreatePrimitiveOpenKnotVector(connect_order);
        pNewFESpace->SetKnotVector(TDim-1, new_knot_vector);
        pNewFESpace->SetInfo(TDim-1, connect_order+1, connect_order);

        // create the new patch
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(-1, pNewFESpace));

        // create the new control point grid
        typedef typename Patch<TDim>::ControlPointType ControlPointType;
        typename StructuredControlGrid<TDim-1, ControlPointType>::Pointer pControlPointGrid1
            = boost::dynamic_pointer_cast<StructuredControlGrid<TDim-1, ControlPointType> >(pPatch1->pControlPointGridFunction()->pControlGrid());
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid2 = pPatch2->pControlPointGridFunction()->pControlGrid();

        //// make a size check, it is not necessary anyway
        if (pControlPointGrid1->size() != pControlPointGrid2->size())
            KRATOS_THROW_ERROR(std::logic_error, "The size of two control point grid are not the same", "")

        // assign data to the new control point grid
        std::vector<std::size_t> new_sizes(TDim);
        for (std::size_t dim = 0; dim < TDim-1; ++dim)
            new_sizes[dim] = pControlPointGrid1->Size(dim);
        new_sizes[TDim-1] = 2;
        typename StructuredControlGrid<TDim, ControlPointType>::Pointer pNewControlPointGrid = StructuredControlGrid<TDim, ControlPointType>::Create(new_sizes);
        for (std::size_t i = 0; i < pControlPointGrid1->size(); ++i)
        {
            pNewControlPointGrid->SetData(i, pControlPointGrid1->GetData(i));
            pNewControlPointGrid->SetData(i + pControlPointGrid1->size(), pControlPointGrid2->GetData(i));
        }
        pNewControlPointGrid->SetName(pControlPointGrid1->Name());

        // assign new control point grid to new patch
        pNewPatch->CreateControlPointGridFunction(pNewControlPointGrid);

        // TODO create other grid function data

        return pNewPatch;
    }

    static int GetDimensionOfGeo(const std::string& fn)
    {
        std::ifstream infile(fn.c_str());
        if(!infile)
            KRATOS_THROW_ERROR(std::logic_error, "Error open file", fn)

        std::string line;
        std::vector<std::string> words;
        int read_mode = _READ_PATCH;
        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == _READ_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                    }

                    // read info
                    int Dim = atoi(words[0].c_str());
                    return Dim;
                }
            }
        }

        return 0;
    }

    /// Create the B-Splines patch from geo file
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchFromGeo(const std::string& fn)
    {
        std::ifstream infile(fn.c_str());
        if(!infile)
            KRATOS_THROW_ERROR(std::logic_error, "Error open file", fn)

        std::string line;
        std::vector<std::string> words;
        int read_mode = _READ_PATCH;
        int npatches, dim_index = 0;
        std::vector<std::size_t> orders;
        std::vector<std::size_t> numbers;
        std::vector<std::vector<double> > knots(3);
        std::vector<std::vector<double> > wcoords(3);
        std::vector<double> weights;
        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == _READ_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                    }

                    // read info
                    int Dim = atoi(words[0].c_str());
                    if (Dim != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The input dimension is invalid", "")
                    npatches = atoi(words[1].c_str());
                    if(npatches > 1)
                    {
                        KRATOS_WATCH(line)
                        KRATOS_WATCH(words[0])
                        KRATOS_WATCH(words[1])
                        KRATOS_THROW_ERROR(std::logic_error, "At present, the number of patches > 1 is not supported, npatches =", npatches)
                    }
                    read_mode = _READ_ORDER;
                    continue;
                }

                if(read_mode == _READ_ORDER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Order section must contained number of information equal to dimension, current number of information =", words.size())

                    // read info
                    for(std::size_t i = 0; i < TDim; ++i)
                        orders.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_NUMBER;
                    continue;
                }

                if(read_mode == _READ_NUMBER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Number section must contained number of information equal to dimension, current number of information =", words.size())

                    for(std::size_t i = 0; i < TDim; ++i)
                        numbers.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_KNOTS;
                    continue;
                }

                if(read_mode == _READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers[dim_index] + orders[dim_index] + 1;
                    if(words.size() != knot_len)
                        KRATOS_THROW_ERROR(std::logic_error, "The Knots section must contained number of information equal to n+p+1, current number of information =", words.size())

                    for(std::size_t i = 0; i < knot_len; ++i)
                    {
                        double k = atof(words[i].c_str());
                        knots[dim_index].push_back(k);
                    }

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = _READ_COORDINATES;
                    }
                    continue;
                }

                if(read_mode == _READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Coordinates section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        wcoords[dim_index].push_back(atof(words[i].c_str()));

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = _READ_WEIGHTS;
                    }
                    continue;
                }

                if(read_mode == _READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Weights section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        weights.push_back(atof(words[i].c_str()));

                    read_mode = _NO_READ;
                    continue;
                }
            }
        }

        // close the file
        infile.close();

        // create the FESpace
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
        for (int dim = 0; dim < TDim; ++dim)
        {
            pNewFESpace->SetKnotVector(dim, knots[dim]);
            pNewFESpace->SetInfo(dim, numbers[dim], orders[dim]);
        }

        // reset function indices and enumerate it first time to give each function in the FESpace a different id
        pNewFESpace->ResetFunctionIndices();
        std::size_t start = 0;
        pNewFESpace->Enumerate(start);

        // create new patch
        typename Patch<TDim>::Pointer pNewPatch = Patch<TDim>::Create(0, pNewFESpace);

        // create control grid and assign to new patch
        typedef ControlPoint<double> ControlPointType;
        typename StructuredControlGrid<TDim, ControlPointType>::Pointer pControlPointGrid = StructuredControlGrid<TDim, ControlPointType>::Create(numbers);
        std::size_t total_number = 1;
        for (int dim = 0; dim < TDim; ++dim)
            total_number *= numbers[dim];

        for (std::size_t i = 0; i < total_number; ++i)
        {
            ControlPointType c;
            if (TDim == 2)
                c.SetCoordinates(wcoords[0][i]/weights[i], wcoords[1][i]/weights[i], 0.0, weights[i]);
            else if (TDim == 3)
                c.SetCoordinates(wcoords[0][i]/weights[i], wcoords[1][i]/weights[i], wcoords[2][i]/weights[i], weights[i]);
            pControlPointGrid->SetData(i, c);
        }

        pControlPointGrid->SetName("CONTROL_POINT");
        pNewPatch->CreateControlPointGridFunction(pControlPointGrid);

        std::cout << __FUNCTION__ << ": Read NURBS from " << fn << " completed" << std::endl;
        return pNewPatch;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplinesPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED defined


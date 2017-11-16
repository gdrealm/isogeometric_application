//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/grid_function.h"
#include "custom_utilities/nurbs/patch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class MultiPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiPatchUtility() {}

    /// Destructor
    virtual ~MultiPatchUtility() {}

    /// Create new patch and wrap it with pointer
    template<int TDim>
    typename Patch<TDim>::Pointer CreatePatchPointer(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace) const
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace));
    }

    /// Export a single NURBS patch in 1D to geometry file that can be read by geo_load routine.
    /// As of 2017, geo_load must be modified to read in the line.
    template<int TDim>
    void ExportGeo(typename Patch<TDim>::Pointer pPatch, const std::string& filename) const
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        if (pPatch->FESpace()->Type() != NURBSFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "does not support non-NURBS patch")

        typename NURBSFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<NURBSFESpace<TDim> >(pPatch->FESpace());

        outfile << "# nurbs mesh v.0.6\n";
        outfile << "#\n";
        outfile << "# NURBS representation for patch " << pPatch->Id() << "\n";
        outfile << "#\n";

        outfile << "#dim\n";
        outfile << TDim << " 1\n";
        outfile << "#p\n";
        for (std::size_t dim = 0; dim < TDim; ++dim)
            outfile << " " << pFESpace->Order(dim);
        outfile << "\n";
        outfile << "#n\n";
        for (std::size_t dim = 0; dim < TDim; ++dim)
            outfile << " " << pFESpace->Number(dim);
        outfile << "\n";

        outfile << "#knots\n";
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            for (std::size_t i = 0; i < pFESpace->KnotVector(dim).size(); ++i)
                outfile << " " << pFESpace->KnotVector(dim)[i];
            outfile << std::endl;
        }

        this->WriteControlPoints<TDim>(outfile, pPatch);

        outfile << std::endl;

        outfile.close();

        std::cout << pPatch->Type() << " " << pPatch->Id() << " is exported to " << filename << " successfully" << std::endl;
    }

    /// Export the multipatch to Glvis for visualization
    template<int TDim>
    void ExportGlvis(typename MultiPatch<TDim>::Pointer pPatch, const std::string& filename) const
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        outfile << "MFEM NURBS mesh v1.0\n\n";

        outfile << "#\n";
        outfile << "# MFEM Geometry Types (see mesh/geom.hpp):\n";
        outfile << "#\n";
        outfile << "# SEGMENT     = 1\n";
        outfile << "# SQUARE      = 3\n";
        outfile << "# CUBE        = 5\n";
        outfile << "#\n\n";

        outfile << "dimension\n" << TDim << "\n\n";

        std::size_t nvertices;
        std::vector<std::vector<std::size_t> > elements;
        std::vector<std::vector<std::size_t> > boundary;
        std::vector<std::size_t> boundary_attr;
        std::vector<std::size_t> edges_knotv;
        // pPatch->GenerateCornerTopology(nvertices, elements, boundary, boundary_attr, edges, edges_knotv);

        outfile << "patches\n\n";
        for (typename MultiPatch<TDim>::PatchContainerType::iterator it = pPatch->begin(); it != pPatch->end(); ++it)
        {
            outfile << "# patch " << it->Id() << "\n\n";
            outfile << "knotvectors\n" << TDim << "\n";

            if (it->FESpace()->Type() != NURBSFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "does not support non-NURBS patch")

            typename NURBSFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<NURBSFESpace<TDim> >(it->FESpace());
            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                outfile << pFESpace->Order(dim) << " " << pFESpace->Number(dim);
                for (std::size_t i = 0; i < pFESpace->KnotVector(dim).size(); ++i)
                    outfile << " " << pFESpace->KnotVector(dim)[i];
                outfile << "\n";
            }
            outfile << "\n";

            outfile << "dimension\n" << TDim << "\n\n";

            typename ControlGrid<ControlPoint<double> >::Pointer pControlGrid = it->ControlPointGridFunction()->ControlGrid();
            outfile << "controlpoints\n";
            for (std::size_t i = 0; i < pControlGrid->size(); ++i)
            {
                for (std::size_t dim = 0; dim < TDim; ++dim)
                    outfile << " " << (*pControlGrid)[i][dim];
                outfile << " " << (*pControlGrid)[i][3] << "\n";
            }
            outfile << "\n";
        }

        outfile << std::endl;

        outfile.close();

        std::cout <<" Multipatch is exported to " << filename << " successfully" << std::endl;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    template<int TDim>
    void WriteControlPoints(std::ostream& rOStream, typename Patch<TDim>::Pointer pPatch) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "WriteControlPoints is not implemented for dimension", TDim)
    }

};

template<>
void MultiPatchUtility::WriteControlPoints<1>(std::ostream& rOStream, Patch<1>::Pointer pPatch) const
{
    typedef Patch<1>::ControlPointType ControlPointType;
    typename RegularControlGrid<1, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const RegularControlGrid<1, ControlPointType> >(pPatch->ControlPointGridFunction()->ControlGrid());

    rOStream << "#u\n";
    for (std::size_t dim = 0; dim < 3; ++dim)
    {
        for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
            rOStream << "\t" << pControlPointGrid->GetValue(nu)(dim);
        rOStream << std::endl;
    }
    for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
        rOStream << "\t" << pControlPointGrid->GetValue(nu)(3);
}

template<>
void MultiPatchUtility::WriteControlPoints<2>(std::ostream& rOStream, Patch<2>::Pointer pPatch) const
{
    typedef Patch<2>::ControlPointType ControlPointType;
    typename RegularControlGrid<2, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const RegularControlGrid<2, ControlPointType> >(pPatch->ControlPointGridFunction()->ControlGrid());

    rOStream << "#u v\n";
    for (std::size_t dim = 0; dim < 2; ++dim)
    {
        for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
            for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
                rOStream << "\t" << pControlPointGrid->GetValue(nu, nv)(dim);
        rOStream << std::endl;
    }
    for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
        for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
            rOStream << "\t" << pControlPointGrid->GetValue(nu, nv)(3);
}

template<>
void MultiPatchUtility::WriteControlPoints<3>(std::ostream& rOStream, Patch<3>::Pointer pPatch) const
{
    typedef Patch<3>::ControlPointType ControlPointType;
    typename RegularControlGrid<3, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const RegularControlGrid<3, ControlPointType> >(pPatch->ControlPointGridFunction()->ControlGrid());

    rOStream << "#u v w\n";
    for (std::size_t dim = 0; dim < 3; ++dim)
    {
        for (std::size_t nw = 0; nw < pControlPointGrid->Size(2); ++nw)
            for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
                for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
                    rOStream << "\t" << pControlPointGrid->GetValue(nu, nv, nw)(dim);
        rOStream << std::endl;
    }
    for (std::size_t nw = 0; nw < pControlPointGrid->Size(2); ++nw)
        for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
            for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
                rOStream << "\t" << pControlPointGrid->GetValue(nu, nv, nw)(3);
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED defined


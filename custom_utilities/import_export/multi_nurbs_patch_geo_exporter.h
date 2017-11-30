//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_EXPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

class MultiNURBSPatchGeoExporterHelper
{
public:
    template<int TDim>
    static void WriteGeoControlPoints(std::ostream& rOStream, typename Patch<TDim>::Pointer pPatch)
    {
        KRATOS_THROW_ERROR(std::logic_error, "WriteGeoControlPoints is not implemented for dimension", TDim)
    }
};

/**
Export NURBS patch/multipatch to Geo to visualize with NURBS toolbox by M. Spink
 */
template<int TDim>
class MultiNURBSPatchGeoExporterWriter : public MultiPatchExporter<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGeoExporterWriter);

    /// Type definition
    typedef MultiPatchExporter<TDim> BaseType;
    typedef typename BaseType::knot_container_t knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiNURBSPatchGeoExporterWriter() : BaseType() {}

    /// Destructor
    virtual ~MultiNURBSPatchGeoExporterWriter() {}

    /// Export a single patch
    virtual void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename) const
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "does not support non-NURBS patch")

        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());

        outfile << "# nurbs mesh v.0.6\n";
        outfile << "#\n";
        outfile << "# BSplines representation for patch " << pPatch->Id() << "\n";
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

        MultiNURBSPatchGeoExporterHelper::WriteGeoControlPoints<TDim>(outfile, pPatch);

        outfile << std::endl;

        outfile.close();

        std::cout << pPatch->Type() << " " << pPatch->Id() << " is exported to " << filename << " successfully" << std::endl;
    }

    /// Export a single patch
    virtual void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename) const
    {
        BaseType::Export(pMultiPatch, filename);
    }
}; // end class MultiNURBSPatchGeoExporterWriter


/// Function instantiation
template<>
void MultiNURBSPatchGeoExporterHelper::WriteGeoControlPoints<1>(std::ostream& rOStream, Patch<1>::Pointer pPatch)
{
    typedef Patch<1>::ControlPointType ControlPointType;
    typename StructuredControlGrid<1, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<1, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());

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


/// Function instantiation
template<>
void MultiNURBSPatchGeoExporterHelper::WriteGeoControlPoints<2>(std::ostream& rOStream, Patch<2>::Pointer pPatch)
{
    typedef Patch<2>::ControlPointType ControlPointType;
    typename StructuredControlGrid<2, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<2, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());

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


/// Function instantiation
template<>
void MultiNURBSPatchGeoExporterHelper::WriteGeoControlPoints<3>(std::ostream& rOStream, Patch<3>::Pointer pPatch)
{
    typedef Patch<3>::ControlPointType ControlPointType;
    typename StructuredControlGrid<3, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<3, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());

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


class MultiNURBSPatchGeoExporter
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGeoExporter);

    template<int TDim>
    static void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename)
    {
        MultiNURBSPatchGeoExporterWriter<TDim> dummy;
        dummy.Export(pPatch, filename);
    }

    template<int TDim>
    static void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename)
    {
        MultiNURBSPatchGeoExporterWriter<TDim> dummy;
        dummy.Export(pMultiPatch, filename);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiNURBSPatchGeoExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiNURBSPatchGeoExporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}


} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_EXPORTER_H_INCLUDED defined


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/grid_function.h"
#include "custom_utilities/nurbs/nurbs_patch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class NURBSPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSPatchUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSPatchUtility() {}

    /// Destructor
    virtual ~NURBSPatchUtility() {}

    /// Export a single NURBS patch in 1D to geometry file that can be read by geo_load routine.
    /// As of 2017, geo_load must be modified to read in the line.
    static void ExportGeo(NURBSPatch<1>::Pointer pPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        outfile << "# nurbs mesh v.0.6\n";
        outfile << "#\n";
        outfile << "# NURBS representation for patch " << pPatch->Id() << "\n";
        outfile << "#\n";

        outfile << "#dim\n";
        outfile << "1 1\n";
        outfile << "#p\n";
        for (std::size_t dim = 0; dim < 1; ++dim)
            outfile << " " << pPatch->Order(dim);
        outfile << "\n";
        outfile << "#n\n";
        for (std::size_t dim = 0; dim < 1; ++dim)
            outfile << " " << pPatch->Number(dim);
        outfile << "\n";

        outfile << "#knots\n";
        for (std::size_t dim = 0; dim < 1; ++dim)
        {
            for (std::size_t i = 0; i < pPatch->KnotVector(dim).size(); ++i)
                outfile << " " << pPatch->KnotVector(dim)[i];
            outfile << std::endl;
        }

        outfile << "#u v\n";
        typename GridFunction<1, NURBSPatch<1>::ControlPointType>::ConstPointer pControlPointGrid = pPatch->ControlPointGrid();
        for (std::size_t dim = 0; dim < 3; ++dim)
        {
            for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
                outfile << "\t" << pControlPointGrid->GetValue(nu)(dim);
            outfile << std::endl;
        }
        for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
            outfile << "\t" << pControlPointGrid->GetValue(nu)(3);
        outfile << std::endl;

        outfile.close();
    }

    /// Export a single NURBS patch in 2D to geometry file that can be read by geo_load routine.
    static void ExportGeo(NURBSPatch<2>::Pointer pPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        outfile << "# nurbs mesh v.0.6\n";
        outfile << "#\n";
        outfile << "# NURBS representation for patch " << pPatch->Id() << "\n";
        outfile << "#\n";

        outfile << "#dim\n";
        outfile << "2 1\n";
        outfile << "#p\n";
        for (std::size_t dim = 0; dim < 2; ++dim)
            outfile << " " << pPatch->Order(dim);
        outfile << "\n";
        outfile << "#n\n";
        for (std::size_t dim = 0; dim < 2; ++dim)
            outfile << " " << pPatch->Number(dim);
        outfile << "\n";

        outfile << "#knots\n";
        for (std::size_t dim = 0; dim < 2; ++dim)
        {
            for (std::size_t i = 0; i < pPatch->KnotVector(dim).size(); ++i)
                outfile << " " << pPatch->KnotVector(dim)[i];
            outfile << std::endl;
        }

        outfile << "#u v\n";
        typename GridFunction<2, NURBSPatch<2>::ControlPointType>::ConstPointer pControlPointGrid = pPatch->ControlPointGrid();
        for (std::size_t dim = 0; dim < 2; ++dim)
        {
            for (std::size_t nv = 0; nv < pPatch->Number(1); ++nv)
                for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
                    outfile << "\t" << pControlPointGrid->GetValue(nu, nv)(dim);
            outfile << std::endl;
        }
        for (std::size_t nv = 0; nv < pPatch->Number(1); ++nv)
            for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
                outfile << "\t" << pControlPointGrid->GetValue(nu, nv)(3);
        outfile << std::endl;

        outfile.close();
    }

    /// Export a single NURBS patch in 3D to geometry file that can be read by geo_load routine.
    static void ExportGeo(NURBSPatch<3>::Pointer pPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        outfile << "# nurbs mesh v.0.6\n";
        outfile << "#\n";
        outfile << "# NURBS representation for patch " << pPatch->Id() << "\n";
        outfile << "#\n";

        outfile << "#dim\n";
        outfile << "3 1\n";
        outfile << "#p\n";
        for (std::size_t dim = 0; dim < 3; ++dim)
            outfile << " " << pPatch->Order(dim);
        outfile << "\n";
        outfile << "#n\n";
        for (std::size_t dim = 0; dim < 3; ++dim)
            outfile << " " << pPatch->Number(dim);
        outfile << "\n";

        outfile << "#knots\n";
        for (std::size_t dim = 0; dim < 3; ++dim)
        {
            for (std::size_t i = 0; i < pPatch->KnotVector(dim).size(); ++i)
                outfile << " " << pPatch->KnotVector(dim)[i];
            outfile << std::endl;
        }

        outfile << "#u v w\n";
        typename GridFunction<3, NURBSPatch<3>::ControlPointType>::ConstPointer pControlPointGrid = pPatch->ControlPointGrid();
        for (std::size_t dim = 0; dim < 3; ++dim)
        {
            for (std::size_t nw = 0; nw < pPatch->Number(2); ++nw)
                for (std::size_t nv = 0; nv < pPatch->Number(1); ++nv)
                    for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
                        outfile << "\t" << pControlPointGrid->GetValue(nu, nv, nw)(dim);
            outfile << std::endl;
        }
        for (std::size_t nw = 0; nw < pPatch->Number(2); ++nw)
            for (std::size_t nv = 0; nv < pPatch->Number(1); ++nv)
                for (std::size_t nu = 0; nu < pPatch->Number(0); ++nu)
                    outfile << "\t" << pControlPointGrid->GetValue(nu, nv, nw)(3);
        outfile << std::endl;

        outfile.close();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NURBSPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_UTILITY_H_INCLUDED defined


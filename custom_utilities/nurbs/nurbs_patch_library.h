//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_LIBRARY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_LIBRARY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/patch.h"
#include "custom_utilities/nurbs/nurbs_patch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class NURBSPatchLibrary
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSPatchLibrary);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSPatchLibrary() {}

    /// Destructor
    virtual ~NURBSPatchLibrary() {}

    /// Generate regular NURBS patch. For 2D, it's rectangle and for 3D it's a cube.
    /// The knot vector only contains 0 and 1, i.e [0 0 0 ... 1 1 1].
    /// All the weights are 1.
    template<std::size_t TDim>
    static typename Patch<TDim>::Pointer CreateRegularPatch(const std::vector<std::size_t>& Orders)
    {
        typename NURBSPatch<TDim>::Pointer pPatch = typename NURBSPatch<TDim>::Pointer(new NURBSPatch<TDim>());

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            knot_container_t knot_vector;
            for (std::size_t i = 0; i < Orders[dim]+1; ++i)
                knot_vector.pCreateKnot(0.0);
            for (std::size_t i = 0; i < Orders[dim]+1; ++i)
                knot_vector.pCreateKnot(1.0);
            pPatch->SetKnotVector(dim, knot_vector);
        }

        return pPatch;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NURBSPatchLibrary";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSPatchLibrary& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_LIBRARY_H_INCLUDED defined


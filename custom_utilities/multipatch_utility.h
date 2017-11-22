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
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"

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

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

}; // end class MultiPatchUtility


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


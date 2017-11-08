//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_REFINEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_REFINEMENT_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/grid_function.h"
#include "custom_utilities/nurbs/patch.h"

namespace Kratos
{

/**
This class represents a single NURBS patch in parametric coordinates.
 */
class NURBSRefinementUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSRefinementUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSRefinementUtility() {}

    /// Destructor
    virtual ~NURBSRefinementUtility() {}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NURBSRefinementUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSRefinementUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_REFINEMENT_UTILITY_H_INCLUDED defined


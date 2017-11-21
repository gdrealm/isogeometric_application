#include "custom_utilities/hierarchical_bsplines/hb_cell.h"
#include "custom_utilities/hierarchical_bsplines/hb_basis_function.h"

namespace Kratos
{

    void HBCell::PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    void HBCell::PrintData(std::ostream& rOStream) const
    {
        rOStream << ", supporting basis functions: (";
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            rOStream << " " << (*it)->Id();
        rOStream << ")";
        BaseType::PrintData(rOStream);
    }

}


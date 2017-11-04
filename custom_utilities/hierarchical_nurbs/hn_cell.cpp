#include "custom_utilities/hierarchical_nurbs/hn_cell.h"
#include "custom_utilities/hierarchical_nurbs/hn_basis_function.h"

namespace Kratos
{

    void HnCell::PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    void HnCell::PrintData(std::ostream& rOStream) const
    {
        rOStream << ", supporting basis functions: (";
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            rOStream << " " << (*it)->Id();
        rOStream << ")";
        BaseType::PrintData(rOStream);
    }

}


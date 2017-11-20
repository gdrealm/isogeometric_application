//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_INDEXING_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_INDEXING_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class provides sub-routines to index the NURBS basis functions in 1D, 2D, 3D.
The base index is assumed to be 1.
 */
class NURBSIndexingUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSIndexingUtility);

    /// Type definition

    /// Default constructor
    NURBSIndexingUtility() {}

    /// Destructor
    virtual ~NURBSIndexingUtility() {}

    /// Compute the NURBS index
    template<int TDim>
    static std::size_t Index(const std::vector<std::size_t>& index, const std::vector<std::size_t>& size)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Index is not implemented for dimension", TDim)
    }

    /// Compute the NURBS index in 0D.
    static std::size_t Index0D(const std::size_t& i)
    {
        return 0;
    }

    /// Compute the NURBS index in 1D.
    static std::size_t Index1D(const std::size_t& i, const std::size_t& n1)
    {
        return i-1;
    }

    /// Compute the NURBS index in 2D
    static std::size_t Index2D(const std::size_t& i, const std::size_t& j,
        const std::size_t& n1, const std::size_t& n2)
    {
        return (j-1)*n1 + (i-1);
    }

    /// Compute the NURBS index in 3D
    static std::size_t Index3D(const std::size_t& i, const std::size_t& j, const std::size_t& k, const std::size_t& n1, const std::size_t& n2, const std::size_t& n3)
    {
        return ((k-1)*n2 + (j-1))*n1 + (i-1);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NURBSIndexingUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

template<>
std::size_t NURBSIndexingUtility::Index<1>(const std::vector<std::size_t>& index, const std::vector<std::size_t>& size)
{
    return Index1D(index[0], size[0]);
}

template<>
std::size_t NURBSIndexingUtility::Index<2>(const std::vector<std::size_t>& index, const std::vector<std::size_t>& size)
{
    return Index2D(index[0], index[1], size[0], size[1]);
}

template<>
std::size_t NURBSIndexingUtility::Index<3>(const std::vector<std::size_t>& index, const std::vector<std::size_t>& size)
{
    return Index3D(index[0], index[1], index[2], size[0], size[1], size[2]);
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSIndexingUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_INDEXING_UTILITY_H_INCLUDED defined


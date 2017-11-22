//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class provides sub-routines to index the BSplines basis function in 1D, 2D, 3D.
The base index is assumed to be 1.
 */
class BSplinesIndexingUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesIndexingUtility);

    /// Type definition

    /// Default constructor
    BSplinesIndexingUtility() {}

    /// Destructor
    virtual ~BSplinesIndexingUtility() {}

    /// Compute the BSplines index
    template<int TDim>
    static std::size_t Index(const std::vector<std::size_t>& index, const std::vector<std::size_t>& size)
    {
        if (TDim == 1)
            return Index1D(index[0], size[0]);
        else if (TDim == 2)
            return Index2D(index[0], index[1], size[0], size[1]);
        else if (TDim == 3)
            return Index3D(index[0], index[1], index[2], size[0], size[1], size[2]);
        else
            KRATOS_THROW_ERROR(std::logic_error, "Index is not implemented for dimension", TDim)
    }

    /// Compute the BSplines index in 0D.
    static std::size_t Index0D(const std::size_t& i)
    {
        return 0;
    }

    /// Compute the BSplines index in 1D.
    static std::size_t Index1D(const std::size_t& i, const std::size_t& n1)
    {
        return i-1;
    }

    /// Compute the BSplines index in 2D
    static std::size_t Index2D(const std::size_t& i, const std::size_t& j,
        const std::size_t& n1, const std::size_t& n2)
    {
        return (j-1)*n1 + (i-1);
    }

    /// Compute the BSplines index array in 2D
    static std::vector<std::size_t> IndexArray2D(const std::size_t& i,
        const std::size_t& n1, const std::size_t& n2)
    {
        std::vector<std::size_t> loc(2);
        loc[0] = (i-1) % n1 + 1;
        loc[1] = (i-1) / n1 + 1;
        return loc;
    }

    /// Compute the BSplines index in 3D
    static std::size_t Index3D(const std::size_t& i, const std::size_t& j, const std::size_t& k,
        const std::size_t& n1, const std::size_t& n2, const std::size_t& n3)
    {
        return ((k-1)*n2 + (j-1))*n1 + (i-1);
    }

    /// Compute the BSplines index array in 3D
    static std::vector<std::size_t> IndexArray3D(const std::size_t& i,
        const std::size_t& n1, const std::size_t& n2, const std::size_t& n3)
    {
        std::vector<std::size_t> loc(3);
        loc[0] = (i-1) % n1 + 1;
        loc[1] = ((i-1) / n1) % n2 + 1;
        loc[2] = ((std::size_t)((i-1) / n1)) / n2 + 1;
        return loc;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplinesIndexingUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesIndexingUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED defined


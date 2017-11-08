//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class represents a NURBS GridFunction in parametric coordinates.
 */
template<std::size_t TDim, typename TDataType>
class GridFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    /// Type definition

    /// Default constructor
    GridFunction() : mName("UNKNOWN") {}

    /// Constructor with name
    GridFunction(const std::string& Name) : mName(Name) {}

    /// Destructor
    virtual ~GridFunction() {}

    void SetName(const std::string& Name) {mName = Name;}
    const std::string& Name() const {return mName;}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::string mName;
};

template<typename TDataType>
class GridFunction<1, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    /// Default constructor
    GridFunction(const std::size_t& n) : mSize(n)
    {
        mData.resize(n);
        std::fill(mData.begin(), mData.end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size() {return mSize;}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i)
    {
        return mData[i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i) const
    {
        return mData[i];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const TDataType& value)
    {
        mData[i] = value;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "[" << mSize << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ": (";
        for (std::size_t i = 0; i < mData.size(); ++i)
            rOStream << " " << mData[i];
        rOStream << ")" << std::endl;
    }

private:
    std::size_t mSize;
    std::vector<TDataType> mData;
};

template<typename TDataType>
class GridFunction<2, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    /// Default constructor
    GridFunction(const std::size_t& m, const std::size_t& n) : mSize{m, n}
    {
        mData.resize(m*n);
        std::fill(mData.begin(), mData.end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of the grid function
    const std::size_t& Size() {return mSize[0]*mSize[1];}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size(const std::size_t& dim) {return mSize[dim];}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i, const std::size_t& j)
    {
        return mData[i*mSize[1] + j];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i, const std::size_t& j) const
    {
        return mData[i*mSize[1] + j];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const std::size_t& j, const TDataType& value)
    {
        mData[i*mSize[1] + j] = value;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "[" << mSize[0] << ", " << mSize[1] << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ": (";
        for (std::size_t i = 0; i < mSize[0]; ++i)
        {
            rOStream << " (";
            for (std::size_t j = 0; j < mSize[1]; ++j)
            {
                rOStream << " " << GetValue(i, j);
            }
            rOStream << ")" << std::endl;
        }
        rOStream << " )" << std::endl;
    }

private:
    std::size_t mSize[2];
    std::vector<TDataType> mData;
};

template<typename TDataType>
class GridFunction<3, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    /// Default constructor
    GridFunction(const std::size_t& m, const std::size_t& n, const std::size_t& p) : mSize{m, n, p}
    {
        mData.resize(m*n*p);
        std::fill(mData.begin(), mData.end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of the grid function
    const std::size_t& Size() {return mSize[0]*mSize[1]*mSize[2];}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size(const std::size_t& dim) {return mSize[dim];}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k)
    {
        return mData[(i*mSize[1] + j)*mSize[2] + k];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k) const
    {
        return mData[(i*mSize[1] + j)*mSize[2] + k];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k, const TDataType& value)
    {
        mData[(i*mSize[1] + j)*mSize[2] + k] = value;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "[" << mSize[0] << ", " << mSize[1] << ", " << mSize[2] << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ": (";
        for (std::size_t i = 0; i < mSize[0]; ++i)
        {
            rOStream << " (";
            for (std::size_t j = 0; j < mSize[1]; ++j)
            {
                rOStream << " (";
                for (std::size_t k = 0; k < mSize[2]; ++k)
                {
                    rOStream << " " << GetValue(i, j, k);
                }
            }
            rOStream << ")" << std::endl;
        }
        rOStream << " )" << std::endl;
    }

private:
    std::size_t mSize[3];
    std::vector<TDataType> mData;
};

/// output stream function
template<std::size_t TDim, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const GridFunction<TDim, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED defined


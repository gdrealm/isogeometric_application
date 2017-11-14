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
    typedef std::vector<TDataType> DataContainerType;

    /// Default constructor
    GridFunction() : mName("UNKNOWN") {}

    /// Constructor with name
    GridFunction(const std::string& Name) : mName(Name) {}

    /// Destructor
    virtual ~GridFunction() {}

    /// Get and Set the name
    void SetName(const std::string& Name) {mName = Name;}
    const std::string Name() const {return mName;}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const GridFunction<TDim, TDataType>& rOther) {}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const typename GridFunction<TDim, TDataType>::Pointer pOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(GridFunction<TDim, TDataType>& rOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename GridFunction<TDim, TDataType>::Pointer pOther) {}

    /// Clone this grid function
    virtual typename GridFunction<TDim, TDataType>::Pointer Clone() {return NULL;}

    /// Get the size of underlying data
    std::size_t Size() const {return mData.size();}

    /// Get the size of underlying data
    std::size_t size() const {return mData.size();}

    /// resize the underlying container
    void resize(const std::size_t& new_size) {mData.resize(new_size);}

    /// Access the underlying data
    const DataContainerType& Data() const {return mData;}

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    void SetData(const std::size_t& i, const TDataType& value) {mData[i] = value;}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

protected:
    /// Access the underlying data
    DataContainerType& Data() {return mData;}

private:
    std::string mName;
    DataContainerType mData;
};

template<typename TDataType>
class GridFunction<1, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    // type definitions
    typedef GridFunction<0, TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;

    /// Constructor with size
    GridFunction(const std::vector<std::size_t>& sizes) : BaseType(), mSize(sizes[0])
    {
        BaseType::Data().resize(sizes[0]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    GridFunction(const std::size_t& n) : BaseType(), mSize(n)
    {
        BaseType::Data().resize(n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size(const std::size_t& dim) const {return mSize;}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i)
    {
        return BaseType::Data()[i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i) const
    {
        return BaseType::Data()[i];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const TDataType& value)
    {
        BaseType::Data()[i] = value;
    }

    // overload operator []
    TDataType& operator[] (const std::size_t& i)
    {
        return BaseType::Data()[i];
    }

    // overload operator []
    const TDataType& operator[] (const std::size_t& i) const
    {
        return BaseType::Data()[i];
    }

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const GridFunction<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
            KRATOS_THROW_ERROR(std::logic_error, "The size of the grid function is incompatible", "")
        for (std::size_t i = 0; i < this->Size(); ++i)
            this->SetValue(i, rOther.GetValue(i));
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const typename GridFunction<1, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const GridFunction<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
            BaseType::Data().resize(rOther.Size());
        for (std::size_t i = 0; i < this->Size(); ++i)
            this->SetValue(i, rOther.GetValue(i));
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename GridFunction<1, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename GridFunction<1, TDataType>::Pointer Clone() const
    {
        typename GridFunction<1, TDataType>::Pointer pNewGridFunction = typename GridFunction<1, TDataType>::Pointer(new GridFunction<1, TDataType>(mSize));
        pNewGridFunction->CopyFrom(*this);
        return pNewGridFunction;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Grid " << BaseType::Name() << "[" << mSize << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Data:\n (";
        for (std::size_t i = 0; i < BaseType::Data().size(); ++i)
            rOStream << " " << BaseType::Data()[i];
        rOStream << ")" << std::endl;
    }

private:
    std::size_t mSize;
};

template<typename TDataType>
class GridFunction<2, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    // type definitions
    typedef GridFunction<0, TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;

    /// Constructor with size
    GridFunction(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    GridFunction(const std::size_t& m, const std::size_t& n) : mSize{m, n}
    {
        BaseType::Data().resize(m*n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size(const std::size_t& dim) const {return mSize[dim];}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i, const std::size_t& j)
    {
        return BaseType::Data()[j*mSize[0] + i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i, const std::size_t& j) const
    {
        return BaseType::Data()[j*mSize[0] + i];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const std::size_t& j, const TDataType& value)
    {
        BaseType::Data()[j*mSize[0] + i] = value;
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const GridFunction<2, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) ) || ( rOther.Size(1) != this->Size(1) ) )
            KRATOS_THROW_ERROR(std::logic_error, "The size of the grid function is incompatible", "")
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                this->SetValue(i, j, rOther.GetValue(i, j));
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const typename GridFunction<2, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const GridFunction<2, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) ) || ( rOther.Size(1) != this->Size(1) ) )
        {
            BaseType::Data().resize(rOther.Size(0)*rOther.Size(1));
        }
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                this->SetValue(i, j, rOther.GetValue(i, j));
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename GridFunction<2, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename GridFunction<2, TDataType>::Pointer Clone() const
    {
        typename GridFunction<2, TDataType>::Pointer pNewGridFunction = typename GridFunction<2, TDataType>::Pointer(new GridFunction<2, TDataType>(mSize[0], mSize[1]));
        pNewGridFunction->CopyFrom(*this);
        return pNewGridFunction;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Grid " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Data:\n (\n";
        for (std::size_t i = 0; i < mSize[0]; ++i)
        {
            rOStream << "  (";
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
};

template<typename TDataType>
class GridFunction<3, TDataType> : public GridFunction<0, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    // type definitions
    typedef GridFunction<0, TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;

    /// Constructor with size
    GridFunction(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1], sizes[2]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]*sizes[2]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    GridFunction(const std::size_t& m, const std::size_t& n, const std::size_t& p) : mSize{m, n, p}
    {
        BaseType::Data().resize(m*n*p);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~GridFunction() {}

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid function is specific dimension
    const std::size_t& Size(const std::size_t& dim) const {return mSize[dim];}

    /// Get the value at specific grid point
    TDataType& GetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k)
    {
        return BaseType::Data()[(k*mSize[1] + j)*mSize[0] + i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k) const
    {
        return BaseType::Data()[(k*mSize[1] + j)*mSize[0] + i];
    }

    /// Set the value at specific grid point
    void SetValue(const std::size_t& i, const std::size_t& j, const std::size_t& k, const TDataType& value)
    {
        BaseType::Data()[(k*mSize[1] + j)*mSize[0] + i] = value;
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const GridFunction<3, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(0) )
          || ( rOther.Size(1) != this->Size(1) )
          || ( rOther.Size(2) != this->Size(2) ) )
            KRATOS_THROW_ERROR(std::logic_error, "The size of the grid function is incompatible", "")
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                for (std::size_t k = 0; k < this->Size(2); ++k)
                    this->SetValue(i, j, k, rOther.GetValue(i, j, k));
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const typename GridFunction<3, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const GridFunction<3, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) )
          || ( rOther.Size(1) != this->Size(1) )
          || ( rOther.Size(2) != this->Size(2) ) )
            BaseType::Data().resize(rOther.Size(0)*rOther.Size(1)*rOther.Size(2));
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                for (std::size_t k = 0; k < this->Size(2); ++k)
                    this->SetValue(i, j, k, rOther.GetValue(i, j, k));
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename GridFunction<3, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename GridFunction<3, TDataType>::Pointer Clone() const
    {
        typename GridFunction<3, TDataType>::Pointer pNewGridFunction = typename GridFunction<3, TDataType>::Pointer(new GridFunction<3, TDataType>(mSize[0], mSize[1], mSize[2]));
        pNewGridFunction->CopyFrom(*this);
        return pNewGridFunction;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Grid " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << ", " << mSize[2] << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Data:\n (";
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

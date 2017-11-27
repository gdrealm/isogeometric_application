//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_grid.h"

namespace Kratos
{

/**
Base class for control value container by a regular grid
*/
template<typename TDataType>
class BaseRegularControlGrid : public ControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BaseRegularControlGrid);

    /// Type definition
    typedef ControlGrid<TDataType> BaseType;
    typedef std::vector<TDataType> DataContainerType;

    /// Default constructor
    BaseRegularControlGrid() : BaseType() {}

    /// Constructor with name
    BaseRegularControlGrid(const std::string& Name) : BaseType(Name) {}

    /// Destructor
    virtual ~BaseRegularControlGrid() {}

    /************************************/
    /********* INHERIT UPSTREAM *********/
    /************************************/

    /// Get the size of underlying data
    virtual std::size_t Size() const {return mData.size();}

    /// Get the size of underlying data
    virtual std::size_t size() const {return mData.size();}

    /// Get the data at specific point
    virtual const TDataType& GetData(const std::size_t& i) const {return mData[i];}

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    virtual void SetData(const std::size_t& i, const TDataType& value) {mData[i] = value;}

    /// overload operator []
    virtual TDataType& operator[] (const std::size_t& i) {return mData[i];}

    /// overload operator []
    virtual const TDataType& operator[] (const std::size_t& i) const {return mData[i];}

    /************************************/
    /****** EXCLUSIVE SUBROUTINES *******/
    /************************************/

    /// resize the underlying container
    void Resize(const std::size_t& new_size) {mData.resize(new_size);}

    /// resize the underlying container
    void resize(const std::size_t& new_size) {mData.resize(new_size);}

    /// Access the underlying data
    DataContainerType& Data() {return mData;}

    /// Access the underlying data
    const DataContainerType& Data() const {return mData;}

    /************************************/
    /******** SUCCEED DOWNSTREAM ********/
    /************************************/

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const ControlGrid<TDataType>& rOther) {}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(ControlGrid<TDataType>& rOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) {}

private:

    DataContainerType mData;
};


/**
Class for control value container by a regular grid
*/
template<int TDim, typename TDataType>
class RegularControlGrid : public BaseRegularControlGrid<TDataType>
{
};


template<typename TDataType>
class RegularControlGrid<1, TDataType> : public BaseRegularControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(RegularControlGrid);

    // type definitions
    typedef BaseRegularControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    RegularControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize(sizes[0])
    {
        BaseType::Data().resize(sizes[0]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    RegularControlGrid(const std::size_t& n) : BaseType(), mSize(n)
    {
        BaseType::Data().resize(n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~RegularControlGrid() {}

    /// Create a new control grid pointer
    static typename RegularControlGrid<1, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename RegularControlGrid<1, TDataType>::Pointer(new RegularControlGrid<1, TDataType>(sizes));
    }

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

    // overload operator ()
    TDataType& operator() (const std::size_t& i) {return BaseType::Data()[i];}

    // overload operator ()
    const TDataType& operator() (const std::size_t& i) const {return BaseType::Data()[i];}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const RegularControlGrid<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
            KRATOS_THROW_ERROR(std::logic_error, "The size of the grid function is incompatible", "")
        for (std::size_t i = 0; i < this->Size(); ++i)
            this->SetValue(i, rOther.GetValue(i));
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const typename RegularControlGrid<1, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const RegularControlGrid<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
            BaseType::Data().resize(rOther.Size());
        for (std::size_t i = 0; i < this->Size(); ++i)
            this->SetValue(i, rOther.GetValue(i));
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename RegularControlGrid<1, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename RegularControlGrid<1, TDataType>::Pointer Clone() const
    {
        typename RegularControlGrid<1, TDataType>::Pointer pNewControlGrid = typename RegularControlGrid<1, TDataType>::Pointer(new RegularControlGrid<1, TDataType>(mSize));
        pNewControlGrid->CopyFrom(*this);
        return pNewControlGrid;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RegularGrid<1> " << BaseType::Name() << "[" << mSize << "]";
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
class RegularControlGrid<2, TDataType> : public BaseRegularControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(RegularControlGrid);

    // type definitions
    typedef BaseRegularControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    RegularControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    RegularControlGrid(const std::size_t& m, const std::size_t& n) : mSize{m, n}
    {
        BaseType::Data().resize(m*n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~RegularControlGrid() {}

    /// Create a new control grid pointer
    static typename RegularControlGrid<2, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename RegularControlGrid<2, TDataType>::Pointer(new RegularControlGrid<2, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(const std::size_t& new_size1, const std::size_t& new_size2)
    {
        resize(new_size1, new_size2);
    }

    /// resize the grid
    void resize(const std::size_t& new_size1, const std::size_t& new_size2)
    {
        mSize[0] = new_size1;
        mSize[1] = new_size2;
        BaseType::Data().resize(new_size1*new_size2);
    }

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

    // overload operator ()
    TDataType& operator() (const std::size_t& i, const std::size_t& j) {return BaseType::Data()[j*mSize[0] + i];}

    // overload operator ()
    const TDataType& operator() (const std::size_t& i, const std::size_t& j) const {return BaseType::Data()[j*mSize[0] + i];}

    /// Copy the data the other grid function
    virtual void CopyFrom(const RegularControlGrid<2, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) ) || ( rOther.Size(1) != this->Size(1) ) )
            KRATOS_THROW_ERROR(std::logic_error, "The size of the grid function is incompatible", "")
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                this->SetValue(i, j, rOther.GetValue(i, j));
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const typename RegularControlGrid<2, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const RegularControlGrid<2, TDataType>& rOther)
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
    virtual void ResizeAndCopyFrom(const typename RegularControlGrid<2, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename RegularControlGrid<2, TDataType>::Pointer Clone() const
    {
        typename RegularControlGrid<2, TDataType>::Pointer pNewControlGrid = typename RegularControlGrid<2, TDataType>::Pointer(new RegularControlGrid<2, TDataType>(mSize[0], mSize[1]));
        pNewControlGrid->CopyFrom(*this);
        return pNewControlGrid;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RegularGrid<2> " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << "]";
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
class RegularControlGrid<3, TDataType> : public BaseRegularControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(RegularControlGrid);

    // type definitions
    typedef BaseRegularControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    RegularControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1], sizes[2]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]*sizes[2]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Constructor with size
    RegularControlGrid(const std::size_t& m, const std::size_t& n, const std::size_t& p) : mSize{m, n, p}
    {
        BaseType::Data().resize(m*n*p);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType(0.0));
    }

    /// Destructor
    virtual ~RegularControlGrid() {}

    /// Create a new control grid pointer
    static typename RegularControlGrid<3, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename RegularControlGrid<3, TDataType>::Pointer(new RegularControlGrid<3, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(const std::size_t& new_size1, const std::size_t& new_size2, const std::size_t& new_size3)
    {
        resize(new_size1, new_size2, new_size3);
    }

    /// resize the grid
    void resize(const std::size_t& new_size1, const std::size_t& new_size2, const std::size_t& new_size3)
    {
        mSize[0] = new_size1;
        mSize[1] = new_size2;
        mSize[2] = new_size3;
        BaseType::Data().resize(new_size1*new_size2*new_size3);
    }

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

    // overload operator ()
    TDataType& operator() (const std::size_t& i, const std::size_t& j, const std::size_t& k)
    {
        return BaseType::Data()[(k*mSize[1] + j)*mSize[0] + i];
    }

    // overload operator ()
    const TDataType& operator() (const std::size_t& i, const std::size_t& j, const std::size_t& k) const
    {
        return BaseType::Data()[(k*mSize[1] + j)*mSize[0] + i];
    }

    /// Copy the data the other grid function
    virtual void CopyFrom(const RegularControlGrid<3, TDataType>& rOther)
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
    virtual void CopyFrom(const typename RegularControlGrid<3, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const RegularControlGrid<3, TDataType>& rOther)
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
    virtual void ResizeAndCopyFrom(const typename RegularControlGrid<3, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Clone this grid function
    virtual typename RegularControlGrid<3, TDataType>::Pointer Clone() const
    {
        typename RegularControlGrid<3, TDataType>::Pointer pNewControlGrid = typename RegularControlGrid<3, TDataType>::Pointer(new RegularControlGrid<3, TDataType>(mSize[0], mSize[1], mSize[2]));
        pNewControlGrid->CopyFrom(*this);
        return pNewControlGrid;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RegularGrid<3> " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << ", " << mSize[2] << "]";
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
template<int TDim, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const RegularControlGrid<TDim, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED defined

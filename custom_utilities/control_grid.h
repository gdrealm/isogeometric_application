//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class is a container to keep the control values.
 */
template<typename TDataType>
class ControlGrid
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlGrid);

    /// Type definition
    typedef TDataType DataType;
    typedef std::vector<TDataType> DataContainerType;

    /// Default constructor
    ControlGrid() : mName("UNKNOWN") {}

    /// Constructor with name
    ControlGrid(const std::string& Name) : mName(Name) {}

    /// Destructor
    virtual ~ControlGrid() {}

    /// Get and Set the name
    void SetName(const std::string& Name) {mName = Name;}
    const std::string Name() const {return mName;}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const ControlGrid<TDataType>& rOther) {}

    /// Copy the data the other grid function. The size of two grid functions must be equal.
    virtual void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(ControlGrid<TDataType>& rOther) {}

    /// Copy the data the other grid function. In the case that the source has different size, the grid function is resized.
    virtual void ResizeAndCopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) {}

    /// Clone this grid function
    virtual typename ControlGrid<TDataType>::Pointer Clone() {return NULL;}

    /// Get the size of underlying data
    std::size_t Size() const {return mData.size();}

    /// Get the size of underlying data
    std::size_t size() const {return mData.size();}

    /// resize the underlying container
    void resize(const std::size_t& new_size) {mData.resize(new_size);}

    /// Access the underlying data
    const DataContainerType& Data() const {return mData;}

    /// Get the data at specific point
    const TDataType& GetData(const std::size_t& i) const {return mData[i];}

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    void SetData(const std::size_t& i, const TDataType& value) {mData[i] = value;}

    // overload operator []
    TDataType& operator[] (const std::size_t& i) {return mData[i];}

    // overload operator []
    const TDataType& operator[] (const std::size_t& i) const {return mData[i];}

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

/// output stream function
template<typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlGrid<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED defined

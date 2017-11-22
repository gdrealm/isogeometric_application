//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/enable_shared_from_this.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"


namespace Kratos
{

enum BoundarySide
{
    _LEFT_   = 0,
    _RIGHT_  = 1,
    _TOP_    = 3,
    _BOTTOM_ = 2,
    _FRONT_  = 4,
    _BACK_   = 5,
    _NUMBER_OF_BOUNDARY_SIDE = 6
};

/**
An FESpace is a collection of shape function defined over the parametric domain. An isogeometric FESpace can be a NURBS FESpace, a hierarchical NURBS FESpace, or a T-Splines FESpace.
 */
template<int TDim>
class FESpace
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Type definition

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace()
    {
    }

    /// Helper to create new BSplinesFESpace pointer
    static FESpace<TDim>::Pointer Create()
    {
        return FESpace<TDim>::Pointer(new FESpace());
    }

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the string representing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the FESpace
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "FESpace" << TDim << "D";
        return ss.str();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the values of the basis function i at point xi
    virtual double GetValue(const std::size_t& i, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the values of the basis functions at point xi
    virtual std::vector<double> GetValue(const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Reset all the dof numbers for each grid function to -1
    void ResetFunctionIndices()
    {
        if (mFunctionIds.size() != this->TotalNumber())
            mFunctionIds.resize(this->TotalNumber());
        std::fill(mFunctionIds.begin(), mFunctionIds.end(), -1);
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    std::size_t& Enumerate(std::size_t& start)
    {
        for (std::size_t i = 0; i < mFunctionIds.size(); ++i)
        {
            if (mFunctionIds[i] == -1)
                mFunctionIds[i] = start++;
        }

        return start;
    }

    /// Access the function indices
    const std::vector<std::size_t>& FunctionIndices() const {return mFunctionIds;}

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<TDim>& rFESpace1, const BoundarySide& side1,
            const FESpace<TDim>& rFESpace2, const BoundarySide& side2) const
    {
        typename FESpace<TDim-1>::Pointer BFESpace1 = rFESpace1.ConstructBoundaryFESpace(side1);
        typename FESpace<TDim-1>::Pointer BFESpace2 = rFESpace1.ConstructBoundaryFESpace(side2);

        return (*BFESpace1) == (*BFESpace2);
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Construct the boundary FESpace based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two FESpacees in terms of its parametric information.
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        return false;
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpace<TDim>& rOther) const
    {
        return this->IsCompatible(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Fast function to get the opposite boundary side
    static BoundarySide OppositeBoundarySide(const BoundarySide& side)
    {
        if (side == _LEFT_) return _RIGHT_;
        else if (side == _RIGHT_) return _LEFT_;
        else if (side == _TOP_) return _BOTTOM_;
        else if (side == _BOTTOM_) return _TOP_;
        else if (side == _FRONT_) return _BACK_;
        else if (side == _BACK_) return _FRONT_;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid boundary side", side)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Overload assignment operator
    FESpace<TDim>& operator=(const FESpace<TDim>& rOther)
    {
        this->mFunctionIds = rOther.mFunctionIds;
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename FESpace<TDim>::Pointer pNewFESpace = typename FESpace<TDim>::Pointer(new FESpace<TDim>());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Add = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Function Indices:";
        for (std::size_t i = 0; i < mFunctionIds.size(); ++i)
            rOStream << " " << mFunctionIds[i];
    }

protected:

    /**
     * data for grid function interpolation
     */
    std::vector<std::size_t> mFunctionIds; // this is to store a unique number of the shape function over the forest of FESpace(s).

private:

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }
};



/**
 * Template specific instantiation for null-D FESpace to terminate the compilation
 */
template<>
class FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace0D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpace<0>& rOther)
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<0>& rFESpace1, const BoundarySide& side1,
            const FESpace<0>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/**
 * Template specific instantiation for -1-D FESpace to terminate the compilation
 */
template<>
class FESpace<-1>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace<-1>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpace<-1>& rOther)
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<-1>& rFESpace1, const BoundarySide& side1,
            const FESpace<-1>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-2>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<-1>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/**
 * Template specific instantiation for -2-D FESpace to terminate the compilation
 */
template<>
class FESpace<-2>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace<-2>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpace<-2>& rOther)
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<-2>& rFESpace1, const BoundarySide& side1,
            const FESpace<-2>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-3>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<-2>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const FESpace<TDim>& rThis)
{
    rOStream << "-------------Begin FESpaceInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End FESpaceInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED defined


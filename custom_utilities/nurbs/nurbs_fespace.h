//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/fespace.h"

namespace Kratos
{

/**
This class represents the FESpace for a single NURBS patch defined over parametric domain.
 */
template<int TDim>
class NURBSFESpace : public FESpace<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSFESpace);

    /// Type definition
    typedef FESpace<TDim> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSFESpace() : BaseType() {}

    /// Destructor
    virtual ~NURBSFESpace() {}

    /// Get the order of the NURBS patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return mOrders[i];}

    /// Get the number of control points of the NURBSFESpace in specific direction
    const std::size_t Number(const std::size_t& i) const {return mNumbers[i];}

    /// Get the number of basis functions defined over the NURBS NURBSFESpace
    virtual const std::size_t TotalNumber() const
    {
        std::size_t Number = 1;
        for (std::size_t i = 0; i < TDim; ++i)
            Number *= mNumbers[i];
        return Number;
    }

    /// Get the string representing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "NURBSFESpace" << TDim << "D";
        return ss.str();
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {
        mpKnotVectors[i] = p_knot_vector;
    }

    /// Create and set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const std::vector<double>& values)
    {
        if (i >= TDim)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
        }
        else
        {
            for (std::size_t j = 0; j < values.size(); ++j)
                mpKnotVectors[i].pCreateKnot(values[j]);
        }
    }

    /// Get the knot vector in i-direction
    const knot_container_t& KnotVector(const std::size_t& i) const {return mpKnotVectors[i];}

    /// Set the NURBS information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {
        mOrders[i] = Order;
        mNumbers[i] = Number;
    }

    /// Validate the NURBSFESpace
    virtual bool Validate() const
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (mpKnotVectors[i].size() != mNumbers[i] + mOrders[i] + 1)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The knot vector is incompatible at dimension", i)
                return false;
            }
        }

        return BaseType::Validate();
    }

    /// Compare between two NURBS patches in terms of parametric information
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const NURBSFESpace<TDim>& rOtherNURBSFESpace = dynamic_cast<const NURBSFESpace<TDim>&>(rOtherFESpace);

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Number(i)) == rOtherNURBSFESpace.Number(i))
                return false;
            if (!(this->Order(i)) == rOtherNURBSFESpace.Order(i))
                return false;
            if (!(this->KnotVector(i) == rOtherNURBSFESpace.KnotVector(i)))
                return false;
        }

        return true;
    }

    /// Construct the boundary patch based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typename NURBSFESpace<TDim-1>::Pointer pBFESpace = typename NURBSFESpace<TDim-1>::Pointer(new NURBSFESpace<TDim-1>());

        if (TDim == 2)
        {
            if (side == _LEFT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(1));
                pBFESpace->SetInfo(0, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _RIGHT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(1));
                pBFESpace->SetInfo(0, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _TOP_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                // TODO set the control points and grid functions
            }
            else if (side == _BOTTOM_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                // TODO set the control points and grid functions
            }
        }
        else if (TDim == 3)
        {
            if (side == _LEFT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(1));
                pBFESpace->SetKnotVector(1, KnotVector(2));
                pBFESpace->SetInfo(0, Number(1), Order(1));
                pBFESpace->SetInfo(1, Number(2), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _RIGHT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(1));
                pBFESpace->SetKnotVector(1, KnotVector(2));
                pBFESpace->SetInfo(0, Number(1), Order(1));
                pBFESpace->SetInfo(1, Number(2), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _TOP_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetKnotVector(1, KnotVector(1));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _BOTTOM_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetKnotVector(1, KnotVector(1));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _FRONT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetKnotVector(1, KnotVector(2));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(1), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _BACK_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0));
                pBFESpace->SetKnotVector(1, KnotVector(2));
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(1), Order(2));
                // TODO set the control points and grid functions
            }
        }

        return pBFESpace;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Add = " << this << ", n = (";
        for (std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << this->Number(i);
        rOStream << "), p = (";
        for (std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << this->Order(i);
        rOStream << ")";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            rOStream << " knot vector " << i << ":";
            for (std::size_t j = 0; j < mpKnotVectors[i].size(); ++j)
                rOStream << " " << mpKnotVectors[i].pKnotAt(j)->Value();
            rOStream << std::endl;
        }
    }

private:

    /**
     * internal data to construct the shape functions on the NURBS
     */
    boost::array<std::size_t, TDim> mOrders;
    boost::array<std::size_t, TDim> mNumbers;
    boost::array<knot_container_t, TDim> mpKnotVectors;
};

/**
 * Template specific instantiation for null-D NURBS patch to terminate the compilation
 */
template<>
class NURBSFESpace<0> : public FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSFESpace);

    /// Type definition
    typedef FESpace<0> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSFESpace() : BaseType() {}

    /// Constructor with id
    NURBSFESpace(const std::size_t& Id) : BaseType(Id) {}

    /// Destructor
    virtual ~NURBSFESpace() {}

    /// Get the order of the NURBS patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return 0;}

    /// Get the number of basis functions defined over the NURBS NURBSFESpace
    virtual const std::size_t Number() const {return 0;}

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "NURBSFESpace0D";
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {}

    /// Set the NURBS information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {}

    /// Validate the NURBSFESpace before using
    virtual bool Validate() const
    {
        return BaseType::Validate();
    }

    /// Compare between two NURBS patches in terms of parametric information
    virtual bool IsCompatible(const FESpace<0>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other FESpace type is not " << Type() << std::endl;
            return false;
        }

        return true;
    }
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSFESpace<TDim>& rThis)
{
    rOStream << "-------------Begin NURBSFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << "-------------End NURBSFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_FESPACE_H_INCLUDED defined


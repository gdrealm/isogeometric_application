//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"

namespace Kratos
{

/**
This class represents the FESpace for a single BSplines patch defined over parametric domain.
 */
template<int TDim>
class BSplinesFESpace : public FESpace<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesFESpace);

    /// Type definition
    typedef FESpace<TDim> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    BSplinesFESpace() : BaseType() {}

    /// Destructor
    virtual ~BSplinesFESpace() {}

    /// Helper to create new BSplinesFESpace pointer
    static BSplinesFESpace<TDim>::Pointer Create()
    {
        return BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace());
    }

    /// Get the order of the BSplines patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return mOrders[i];}

    /// Get the number of control points of the BSplinesFESpace in specific direction
    const std::size_t Number(const std::size_t& i) const {return mNumbers[i];}

    /// Get the number of basis functions defined over the BSplines BSplinesFESpace
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
        ss << "BSplinesFESpace" << TDim << "D";
        return ss.str();
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {
        mKnotVectors[i] = p_knot_vector;
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
                mKnotVectors[i].pCreateKnot(values[j]);
        }
    }

    /// Get the knot vector in i-direction
    const knot_container_t& KnotVector(const std::size_t& i) const {return mKnotVectors[i];}

    /// Set the BSplines information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {
        mOrders[i] = Order;
        mNumbers[i] = Number;
    }

    /// Validate the BSplinesFESpace
    virtual bool Validate() const
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (mKnotVectors[i].size() != mNumbers[i] + mOrders[i] + 1)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The knot vector is incompatible at dimension", i)
                return false;
            }
        }

        return BaseType::Validate();
    }

    /// Get the values of the basis function i at point xi
    virtual double GetValue(const std::size_t& i, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "GetValue is not implemented for dimension", TDim)
    }

    /// Get the values of the basis functions at point xi
    virtual std::vector<double> GetValue(const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "GetValue is not implemented for dimension", TDim)
    }

    /// Compare between two BSplines patches in terms of parametric information
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const BSplinesFESpace<TDim>& rOtherBSplinesFESpace = dynamic_cast<const BSplinesFESpace<TDim>&>(rOtherFESpace);

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Number(i)) == rOtherBSplinesFESpace.Number(i))
                return false;
            if (!(this->Order(i)) == rOtherBSplinesFESpace.Order(i))
                return false;
            if (!(this->KnotVector(i) == rOtherBSplinesFESpace.KnotVector(i)))
                return false;
        }

        return true;
    }

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        std::vector<std::size_t> func_indices;

        if (side == _LEFT_)
        {
            if (TDim == 1)
            {
                func_indices.resize(1);
                func_indices[0] = BaseType::mFunctionIds[BSplinesIndexingUtility::Index1D(1, this->Number(0))];
            }
            else if (TDim == 2)
            {
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility::Index1D(j+1, this->Number(1))]
                        = BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(1, j+1, this->Number(0), this->Number(1))];
            }
            else if (TDim == 3)
            {
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility::Index2D(j+1, k+1, this->Number(1), this->Number(2))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(1, j+1, k+1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _RIGHT_)
        {
            if (TDim == 1)
            {
                func_indices.resize(1);
                func_indices[0] = BaseType::mFunctionIds[BSplinesIndexingUtility::Index1D(this->Number(0), this->Number(0))];
            }
            else if (TDim == 2)
            {
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility::Index1D(j+1, this->Number(1))]
                        = BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(this->Number(0), j+1, this->Number(0), this->Number(1))];
            }
            else if (TDim == 3)
            {
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility::Index2D(j+1, k+1, this->Number(1), this->Number(2))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(this->Number(0), j+1, k+1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _FRONT_)
        {
            if (TDim == 2)
            {
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility::Index1D(i+1, this->Number(0))]
                        = BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(i+1, 0, this->Number(0), this->Number(1))];
            }
            else if (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility::Index2D(i+1, k+1, this->Number(0), this->Number(2))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, 0, k+1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BACK_)
        {
            if (TDim == 2)
            {
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility::Index1D(i+1, this->Number(0))]
                        = BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(i+1, this->Number(1), this->Number(0), this->Number(1))];
            }
            else if (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility::Index2D(i+1, k+1, this->Number(0), this->Number(2))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, this->Number(1), k+1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BOTTOM_)
        {
            if (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, j+1, 0, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _TOP_)
        {
            if (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1))]
                            = BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, j+1, this->Number(2), this->Number(0), this->Number(1), this->Number(2))];
            }
        }

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        if (side == _LEFT_)
        {
            if (TDim == 1)
            {
                BaseType::mFunctionIds[BSplinesIndexingUtility::Index1D(1, this->Number(0))] = func_indices[0];
            }
            else if (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(1, j+1, this->Number(0), this->Number(1))]
                        = func_indices[BSplinesIndexingUtility::Index1D(j+1, this->Number(1))];
            }
            else if (TDim == 3)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(1, j+1, k+1, this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(j+1, k+1, this->Number(1), this->Number(2))];
            }
        }
        else if (side == _RIGHT_)
        {
            if (TDim == 1)
            {
                BaseType::mFunctionIds[BSplinesIndexingUtility::Index1D(this->Number(0), this->Number(0))] = func_indices[0];
            }
            else if (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(this->Number(0), j+1, this->Number(0), this->Number(1))]
                        = func_indices[BSplinesIndexingUtility::Index1D(j+1, this->Number(1))];
            }
            else if (TDim == 3)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(this->Number(0), j+1, k+1, this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(j+1, k+1, this->Number(1), this->Number(2))];
            }
        }
        else if (side == _FRONT_)
        {
            if (TDim == 2)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(i+1, 0, this->Number(0), this->Number(1))]
                        = func_indices[BSplinesIndexingUtility::Index1D(i+1, this->Number(0))];
            }
            else if (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, 0, k+1, this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(i+1, k+1, this->Number(0), this->Number(2))];
            }
        }
        else if (side == _BACK_)
        {
            if (TDim == 2)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(i+1, this->Number(1), this->Number(0), this->Number(1))]
                            = func_indices[BSplinesIndexingUtility::Index1D(i+1, this->Number(0))];
            }
            else if (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, this->Number(1), k+1, this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(i+1, k+1, this->Number(0), this->Number(2))];
            }
        }
        else if (side == _BOTTOM_)
        {
            if (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, j+1, 0, this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1))];
            }
        }
        else if (side == _TOP_)
        {
            if (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, j+1, this->Number(2), this->Number(0), this->Number(1), this->Number(2))]
                            = func_indices[BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1))];
            }
        }
    }

    /// Construct the boundary patch based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typename BSplinesFESpace<TDim-1>::Pointer pBFESpace = typename BSplinesFESpace<TDim-1>::Pointer(new BSplinesFESpace<TDim-1>());

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

    /// Overload assignment operator
    BSplinesFESpace<TDim>& operator=(const BSplinesFESpace<TDim>& rOther)
    {
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            this->SetKnotVector(dim, rOther.KnotVector(dim));
            this->SetInfo(dim, rOther.Number(dim), rOther.Order(dim));
        }
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());
        *pNewFESpace = *this;
        return pNewFESpace;
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
            for (std::size_t j = 0; j < mKnotVectors[i].size(); ++j)
                rOStream << " " << mKnotVectors[i].pKnotAt(j)->Value();
            rOStream << std::endl;
        }
        if (BaseType::mFunctionIds.size() == this->TotalNumber())
        {
            rOStream << " Function Indices:";
            if (TDim == 1)
            {
                for (std::size_t i = 0; i < BaseType::mFunctionIds.size(); ++i)
                    rOStream << " " << BaseType::mFunctionIds[i];
            }
            else if (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    for (std::size_t i = 0; i < this->Number(0); ++i)
                        rOStream << " " << BaseType::mFunctionIds[BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1))];
                    rOStream << std::endl;
                }
            }
            else if (TDim == 3)
            {
                for (std::size_t k = 0; k < this->Number(2); ++k)
                {
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                    {
                        for (std::size_t i = 0; i < this->Number(0); ++i)
                            rOStream << " " << BaseType::mFunctionIds[BSplinesIndexingUtility::Index3D(i+1, j+1, k+1, this->Number(0), this->Number(1), this->Number(2))];
                        rOStream << std::endl;
                    }
                    rOStream << std::endl;
                }
            }
        }
    }

private:

    /**
     * internal data to construct the shape functions on the BSplines
     */
    boost::array<std::size_t, TDim> mOrders;
    boost::array<std::size_t, TDim> mNumbers;
    boost::array<knot_container_t, TDim> mKnotVectors;
};

/**
 * Template specific instantiation for null-D BSplines patch to terminate the compilation
 */
template<>
class BSplinesFESpace<0> : public FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesFESpace);

    /// Type definition
    typedef FESpace<0> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    BSplinesFESpace() : BaseType() {}

    /// Destructor
    virtual ~BSplinesFESpace() {}

    /// Get the order of the BSplines patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return 0;}

    /// Get the number of basis functions defined over the BSplines BSplinesFESpace
    virtual const std::size_t Number() const {return 0;}

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "BSplinesFESpace0D";
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {}

    /// Set the BSplines information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {}

    /// Validate the BSplinesFESpace before using
    virtual bool Validate() const
    {
        return BaseType::Validate();
    }

    /// Compare between two BSplines patches in terms of parametric information
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
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesFESpace<TDim>& rThis)
{
    rOStream << "-------------Begin BSplinesFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End BSplinesFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED defined


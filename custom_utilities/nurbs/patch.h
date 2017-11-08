//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/grid_function.h"

namespace Kratos
{

enum BoundarySide
{
    _LEFT_,
    _RIGHT_,
    _TOP_,
    _BOTTOM_,
    _FRONT_,
    _BACK_
};

/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<std::size_t TDim>
class Patch
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef GridFunction<TDim, double> DoubleGridFunctionType;
    typedef std::vector<typename DoubleGridFunctionType::Pointer> DoubleGridFunctionContainterType;
    typedef GridFunction<TDim, array_1d<double, 3> > Array1DGridFunctionType;
    typedef std::vector<typename Array1DGridFunctionType::Pointer> Array1DGridFunctionContainerType;
    typedef GridFunction<TDim, Vector> VectorGridFunctionType;
    typedef std::vector<typename VectorGridFunctionType::Pointer> VectorGridFunctionContainerType;
    typedef std::vector<typename Patch<TDim>::Pointer> NeighborPatchContainerType;

    /// Default constructor
    Patch() : mId(0) {}

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Set the Id of this patch
    void SetId(const std::size_t& Id) {mId = Id;}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t& Number() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t& Order(const std::size_t& i) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the string describing the type of the patch
    virtual const std::string Type() const
    {
        std::stringstream ss;
        ss << "Patch" << TDim << "D";
        return ss.str();
    }

    /// Reset all the dof numbers for each grid function to -1
    void ResetDofs()
    {
        for (std::size_t i = 0; i < mGridDofNumbers.size(); ++i)
            mGridDofNumbers[i] = -1;
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    std::size_t& Enumerate(std::size_t& start)
    {
        for (std::size_t i = 0; i < mGridDofNumbers.size(); ++i)
        {
            if (mGridDofNumbers[i] != -1)
                mGridDofNumbers[i] = start++;
        }

        return start;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Set the control point grid
    void SetControlPointGrid(typename GridFunction<TDim, ControlPointType>::Pointer pControlPointGrid)
    {
        mpControlPointGrid = pControlPointGrid;
    }

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::Pointer GetControlPointGrid() {return mpControlPointGrid;}

    /// Get the control point grid
    const typename GridFunction<TDim, ControlPointType>::Pointer GetControlPointGrid() const {return mpControlPointGrid;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for double variable
    void AddDoubleGridFunction(typename DoubleGridFunctionType::Pointer pGridFunction)
    {
        mpDoubleGridFunctions.push_back(pGridFunction);
    }

    /// Get the underlying double grid functions
    DoubleGridFunctionContainterType& DoubleGridFunctions() {return mpDoubleGridFunctions;}

    /// Get the underlying double grid functions
    const DoubleGridFunctionContainterType& DoubleGridFunctions() const {return mpDoubleGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for array_1d variable
    void AddArray1DGridFunction(typename Array1DGridFunctionType::Pointer pGridFunction)
    {
        mpArray1DGridFunctions.push_back(pGridFunction);
    }

    /// Get the underlying array_1d grid functions
    Array1DGridFunctionContainerType& Array1DGridFunctions() {return mpArray1DGridFunctions;}

    /// Get the underlying array_1d grid functions
    const Array1DGridFunctionContainerType& Array1DGridFunctions() const {return mpArray1DGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for vector variable
    void AddVectorGridFunction(typename VectorGridFunctionType::Pointer pGridFunction)
    {
        mpVectorGridFunctions.push_back(pGridFunction);
    }

    /// Get the underlying vector grid functions
    VectorGridFunctionContainerType& VectorGridFunctions() {return mpVectorGridFunctions;}

    /// Get the underlying vector grid functions
    const VectorGridFunctionContainerType& VectorGridFunctions() const {return mpVectorGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Validate the patch
    virtual bool Validate() const
    {
        if (mpControlPointGrid->Size() != this->Number())
            KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        for (typename DoubleGridFunctionContainterType::const_iterator it = DoubleGridFunctions().begin();
                it != DoubleGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->Number())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The double variable grid is incompatible", (*it)->Name())
                return false;
            }
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions().begin();
                it != Array1DGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->Number())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The array_1d variable grid is incompatible", (*it)->Name())
                return false;
            }
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions().begin();
                it != VectorGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->Number())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The vector variable grid is incompatible", (*it)->Name())
                return false;
            }
        }

        // check the compatibility between patch
        if (pLeft() != NULL)
        {
            if (this->Type() != pLeft()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and left neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _LEFT_, *pLeft(), _RIGHT_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and left neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pRight() != NULL)
        {
            if (this->Type() != pRight()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and right neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _RIGHT_, *pLeft(), _LEFT_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and right neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pTop() != NULL)
        {
            if (this->Type() != pTop()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and right neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _TOP_, *pTop(), _BOTTOM_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and top neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pBottom() != NULL)
        {
            if (this->Type() != pBottom()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and bottom neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _BOTTOM_, *pTop(), _TOP_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and bottom neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pFront() != NULL)
        {
            if (this->Type() != pFront()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and front neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _FRONT_, *pFront(), _BACK_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and front neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pBack() != NULL)
        {
            if (this->Type() != pBack()->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and back neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _BACK_, *pBack(), _FRONT_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and back neighbor is incompatible", "")
                    return false;
                }
            }
        }

        return true;
    }

    /// Check the compatibility between boundaries of two patches
    virtual bool CheckBoundaryCompatibility(const Patch<TDim>& rPatch1, const BoundarySide& side1,
            const Patch<TDim>& rPatch2, const BoundarySide& side2) const
    {
        typename Patch<TDim-1>::Pointer BPatch1 = rPatch1.ConstructBoundaryPatch(side1);
        typename Patch<TDim-1>::Pointer BPatch2 = rPatch1.ConstructBoundaryPatch(side2);

        return (*BPatch1) == (*BPatch2);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Construct the boundary patch based on side
    virtual typename Patch<TDim-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    {
        // TODO
        return NULL;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Set the left neighbor
    void pSetLeft(Patch<TDim>::Pointer pLeft) {mpLeft = pLeft;}

    /// Access the left neighbor
    typename Patch<TDim>::Pointer pLeft() {return mpLeft;}
    const typename Patch<TDim>::Pointer pLeft() const {return mpLeft;}

    /// Set the right neighbor
    void pSetRight(Patch<TDim>::Pointer pRight) {mpRight = pRight;}

    /// Access the right neighbor
    typename Patch<TDim>::Pointer pRight() {return mpRight;}
    const typename Patch<TDim>::Pointer pRight() const {return mpRight;}

    /// Set the top neighbor
    void pSetTop(Patch<TDim>::Pointer pTop) {mpTop = pTop;}

    /// Access the top neighbor
    typename Patch<TDim>::Pointer pTop() {return mpTop;}
    const typename Patch<TDim>::Pointer pTop() const {return mpTop;}

    /// Set the bottom neighbor
    void pSetBottom(Patch<TDim>::Pointer pBottom) {mpBottom = pBottom;}

    /// Access the bottom neighbor
    typename Patch<TDim>::Pointer pBottom() {return mpBottom;}
    const typename Patch<TDim>::Pointer pBottom() const {return mpBottom;}

    /// Set the front neighbor
    void pSetFront(Patch<TDim>::Pointer pFront) {mpFront = pFront;}

    /// Access the front neighbor
    typename Patch<TDim>::Pointer pFront() {return mpFront;}
    const typename Patch<TDim>::Pointer pFront() const {return mpFront;}

    /// Set the back neighbor
    void pSetBack(Patch<TDim>::Pointer pBack) {mpBack = pBack;}

    /// Access the back neighbor
    typename Patch<TDim>::Pointer pBack() {return mpBack;}
    const typename Patch<TDim>::Pointer pBack() const {return mpBack;}

    /// Get all the "real" neighbors of this patch
    NeighborPatchContainerType GetNeighbors() const
    {
        NeighborPatchContainerType pNeighbors;
        if (pLeft() != NULL)
            pNeighbors.push_back(pLeft());
        if (pRight() != NULL)
            pNeighbors.push_back(pRight());
        if (pTop() != NULL)
            pNeighbors.push_back(pTop());
        if (pBottom() != NULL)
            pNeighbors.push_back(pBottom());
        if (pFront() != NULL)
            pNeighbors.push_back(pFront());
        if (pBack() != NULL)
            pNeighbors.push_back(pBack());
        return pNeighbors;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare between two patches in terms of parametric and control points
    virtual bool IsEqual(const Patch<TDim>& rPatch1, const Patch<TDim>& rPatch2) const
    {
        return false;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<TDim>& rOther)
    {
        return this->IsEqual(*this, rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<dim=" << TDim << ">: n = " << this->Number();
        rOStream << ", p = (";
        for (std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << this->Order(i);
        rOStream << ")" << std::endl;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    std::size_t mId;

    /**
     * data for grid function interpolation
     */
    std::vector<std::size_t> mGridDofNumbers; // actually it is not d.o.f; this is to store a unique number to identify the grid function on the forest of patches.

    typename GridFunction<TDim, ControlPointType>::Pointer mpControlPointGrid;
    DoubleGridFunctionContainterType mpDoubleGridFunctions;
    Array1DGridFunctionContainerType mpArray1DGridFunctions;
    VectorGridFunctionContainerType mpVectorGridFunctions;

    /**
     * neighboring data
     */
    typename Patch<TDim>::Pointer mpLeft;
    typename Patch<TDim>::Pointer mpRight;
    typename Patch<TDim>::Pointer mpTop;
    typename Patch<TDim>::Pointer mpBottom;
    typename Patch<TDim>::Pointer mpFront;
    typename Patch<TDim>::Pointer mpBack;
};

/**
 * Template specific instantiation for null-D patch to terminate the compilation
 */
template<>
class Patch<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Default constructor
    Patch() : mId(0) {}

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t& Number() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t& Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the patch
    virtual const std::string Type() const
    {
        return "Patch0D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<0>& rOther)
    {
        return true;
    }

    /// Check the compatibility between boundaries of two patches
    virtual bool CheckBoundaryCompatibility(const Patch<0>& rPatch1, const BoundarySide& side1,
            const Patch<0>& rPatch2, const BoundarySide& side2) const
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const Patch<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED defined


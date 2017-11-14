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
#include <boost/enable_shared_from_this.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/grid_function.h"


#define DEBUG_DESTROY

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


// Declaration
template<std::size_t TDim> class Patch;
template<std::size_t TDim> class MultiPatch;


/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<std::size_t TDim>
class Patch : public boost::enable_shared_from_this<Patch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef GridFunction<0, double> GridFunctionBaseType;
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
    virtual ~Patch()
    {
        #ifdef DEBUG_DESTROY
        std::cout << Type() << ", Id = " << Id() << ", Add = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Set the Id of this patch
    void SetId(const std::size_t& Id) {mId = Id;}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
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
    void SetControlPointGrid(typename GridFunction<TDim, ControlPointType>::Pointer pControlPointGrid){CheckSize(*pControlPointGrid); mpControlPointGrid = pControlPointGrid;}

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::Pointer ControlPointGrid() {return mpControlPointGrid;}

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::ConstPointer ControlPointGrid() const {return mpControlPointGrid;}

    /// Get the control point weights vector
    std::vector<double> GetControlWeights() const
    {
        const typename GridFunction<TDim, ControlPointType>::DataContainerType& GridData = ControlPointGrid()->Data();
        std::vector<double> Weights(GridData.size());
        for (std::size_t i = 0; i < GridData.size(); ++i)
            Weights[i] = GridData[i].W();
        return Weights;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for double variable
    void AddDoubleGridFunction(typename DoubleGridFunctionType::Pointer pGridFunction) {CheckSize(*pGridFunction); mpDoubleGridFunctions.push_back(pGridFunction);}

    /// Get the underlying double grid functions
    DoubleGridFunctionContainterType& DoubleGridFunctions() {return mpDoubleGridFunctions;}

    /// Get the underlying double grid functions
    const DoubleGridFunctionContainterType& DoubleGridFunctions() const {return mpDoubleGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for array_1d variable
    void AddArray1DGridFunction(typename Array1DGridFunctionType::Pointer pGridFunction) {CheckSize(*pGridFunction); mpArray1DGridFunctions.push_back(pGridFunction);}

    /// Get the underlying array_1d grid functions
    Array1DGridFunctionContainerType& Array1DGridFunctions() {return mpArray1DGridFunctions;}

    /// Get the underlying array_1d grid functions
    const Array1DGridFunctionContainerType& Array1DGridFunctions() const {return mpArray1DGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Add the grid function for vector variable
    void AddVectorGridFunction(typename VectorGridFunctionType::Pointer pGridFunction) {CheckSize(*pGridFunction); mpVectorGridFunctions.push_back(pGridFunction);}

    /// Get the underlying vector grid functions
    VectorGridFunctionContainerType& VectorGridFunctions() {return mpVectorGridFunctions;}

    /// Get the underlying vector grid functions
    const VectorGridFunctionContainerType& VectorGridFunctions() const {return mpVectorGridFunctions;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the grid function with name.
    /// One must make sure that there is no repetition in grid function name. Otherwise, the first one will be picked in the order: double grid function, array_1d grid function and vector grid function.
    GridFunctionBaseType::Pointer pGetGridFunction(const std::string& name) const
    {
        for (typename DoubleGridFunctionContainterType::const_iterator it = DoubleGridFunctions().begin();
                it != DoubleGridFunctions().end(); ++it)
        {
            if ((*it)->Name() == name)
                return *it;
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions().begin();
                it != Array1DGridFunctions().end(); ++it)
        {
            if ((*it)->Name() == name)
                return *it;
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions().begin();
                it != VectorGridFunctions().end(); ++it)
        {
            if ((*it)->Name() == name)
                return *it;
        }

        return NULL;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Validate the patch
    virtual bool Validate() const
    {
        if (Id() == 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The patch must have an Id", "")
        }

        if (mpControlPointGrid->Size() != this->TotalNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        for (typename DoubleGridFunctionContainterType::const_iterator it = DoubleGridFunctions().begin();
                it != DoubleGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The double variable grid is incompatible", (*it)->Name())
                return false;
            }
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions().begin();
                it != Array1DGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The array_1d variable grid is incompatible", (*it)->Name())
                return false;
            }
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions().begin();
                it != VectorGridFunctions().end(); ++it)
        {
            if ((*it)->Size() != this->TotalNumber())
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

    /// Synchronize the grid function data with the other patch
    void SynchronizeGridFunction(Patch<TDim>& rOtherPatch)
    {
        // synchronize double grid function from the other patch
        for (typename DoubleGridFunctionContainterType::iterator it = rOtherPatch.DoubleGridFunctions().begin();
                it != rOtherPatch.DoubleGridFunctions().end(); ++it)
        {
            GridFunctionBaseType::Pointer pGridFunction = this->pGetGridFunction((*it)->Name());
            if (pGridFunction != NULL)
            {
                pGridFunction->ResizeAndCopyFrom(*it);
            }
            else
            {
                GridFunctionBaseType::Pointer pNewGridFunction = (*it)->Clone();
                this->DoubleGridFunctions().push_back(pNewGridFunction);
            }
        }

        // synchronize array_1d grid function
        for (typename Array1DGridFunctionContainerType::iterator it = rOtherPatch.Array1DGridFunctions().begin();
                it != rOtherPatch.Array1DGridFunctions().end(); ++it)
        {
            GridFunctionBaseType::Pointer pGridFunction = this->pGetGridFunction((*it)->Name());
            if (pGridFunction != NULL)
            {
                pGridFunction->ResizeAndCopyFrom(*it);
            }
            else
            {
                GridFunctionBaseType::Pointer pNewGridFunction = (*it)->Clone();
                this->Array1DGridFunctions().push_back(pNewGridFunction);
            }
        }

        // synchronize vector grid function
        for (typename VectorGridFunctionContainerType::iterator it = rOtherPatch.VectorGridFunctions().begin();
                it != rOtherPatch.VectorGridFunctions().end(); ++it)
        {
            GridFunctionBaseType::Pointer pGridFunction = this->pGetGridFunction((*it)->Name());
            if (pGridFunction != NULL)
            {
                pGridFunction->ResizeAndCopyFrom(*it);
            }
            else
            {
                GridFunctionBaseType::Pointer pNewGridFunction = (*it)->Clone();
                this->VectorGridFunctions().push_back(pNewGridFunction);
            }
        }

        /* do the same thing with the other patch */
        rOtherPatch.SynchronizeGridFunction(*this);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get/Set the left neighbor
    /// ->shared_from_this() is used to avoid the problem with boost python in assigning the shared_ptr to weak_ptr as reported here: https://stackoverflow.com/questions/8233252/boostpython-and-weak-ptr-stuff-disappearing
    void pSetLeft(Patch<TDim>::Pointer pLeft) {mpLeft = pLeft->shared_from_this();}
    Patch<TDim>& Left() {return *pLeft();}
    const Patch<TDim>& Left() const {return *pLeft();}
    // typename Patch<TDim>::Pointer pLeft() {return mpLeft;}
    // const typename Patch<TDim>::Pointer pLeft() const {return mpLeft;}
    typename Patch<TDim>::Pointer pLeft() {return mpLeft.lock();}
    const typename Patch<TDim>::Pointer pLeft() const {return mpLeft.lock();}

    /// Get/Set the right neighbor
    void pSetRight(Patch<TDim>::Pointer pRight) {mpRight = pRight->shared_from_this();}
    Patch<TDim>& Right() {return *pRight();}
    const Patch<TDim>& Right() const {return *pRight();}
    // typename Patch<TDim>::Pointer& pRefRight() {return mpRight;}
    // const typename Patch<TDim>::Pointer& pRefRight() const {return mpRight;}
    // typename Patch<TDim>::Pointer pRight() {return mpRight;}
    // const typename Patch<TDim>::Pointer pRight() const {return mpRight;}
    typename Patch<TDim>::Pointer pRight() {return mpRight.lock();}
    const typename Patch<TDim>::Pointer pRight() const {return mpRight.lock();}

    /// Get/Set the top neighbor
    void pSetTop(Patch<TDim>::Pointer pTop) {mpTop = pTop->shared_from_this();}
    Patch<TDim>& Top() {return *pTop();}
    const Patch<TDim>& Top() const {return *pTop();}
    // typename Patch<TDim>::Pointer pTop() {return mpTop;}
    // const typename Patch<TDim>::Pointer pTop() const {return mpTop;}
    typename Patch<TDim>::Pointer pTop() {return mpTop.lock();}
    const typename Patch<TDim>::Pointer pTop() const {return mpTop.lock();}

    /// Get/Set the bottom neighbor
    void pSetBottom(Patch<TDim>::Pointer pBottom) {mpBottom = pBottom->shared_from_this();}
    Patch<TDim>& Bottom() {return *pBottom();}
    const Patch<TDim>& Bottom() const {return *pBottom();}
    // typename Patch<TDim>::Pointer pBottom() {return mpBottom;}
    // const typename Patch<TDim>::Pointer pBottom() const {return mpBottom;}
    typename Patch<TDim>::Pointer pBottom() {return mpBottom.lock();}
    const typename Patch<TDim>::Pointer pBottom() const {return mpBottom.lock();}

    /// Get/Set the front neighbor
    void pSetFront(Patch<TDim>::Pointer pFront) {mpFront = pFront->shared_from_this();}
    Patch<TDim>& Front() {return *pFront();}
    const Patch<TDim>& Front() const {return *pFront();}
    // typename Patch<TDim>::Pointer pFront() {return mpFront;}
    // const typename Patch<TDim>::Pointer pFront() const {return mpFront;}
    typename Patch<TDim>::Pointer pFront() {return mpFront.lock();}
    const typename Patch<TDim>::Pointer pFront() const {return mpFront.lock();}

    /// Get/Set the back neighbor
    void pSetBack(Patch<TDim>::Pointer pBack) {mpBack = pBack->shared_from_this();}
    Patch<TDim>& Back() {return *pBack();}
    const Patch<TDim>& Back() const {return *pBack();}
    // typename Patch<TDim>::Pointer pBack() {return mpBack;}
    // const typename Patch<TDim>::Pointer pBack() const {return mpBack;}
    typename Patch<TDim>::Pointer pBack() {return mpBack.lock();}
    const typename Patch<TDim>::Pointer pBack() const {return mpBack.lock();}

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

    /// Get/Set the parent multipatch
    void pSetParentMultiPatch(typename MultiPatch<TDim>::Pointer pPatch) {mpParentMultiPatch = pPatch;}
    MultiPatch<TDim>& ParentMultiPatch() {return *pParentMultiPatch();}
    const MultiPatch<TDim>& ParentMultiPatch() const {return *pParentMultiPatch();}
    typename MultiPatch<TDim>::Pointer pParentMultiPatch() {return mpParentMultiPatch.lock();}
    const typename MultiPatch<TDim>::Pointer pParentMultiPatch() const {return mpParentMultiPatch.lock();}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two patches in terms of its parametric information. The grid function data, including control points, shall not be checked.
    virtual bool IsCompatible(const Patch<TDim>& rOtherPatch) const
    {
        return false;
    }

    /// Compare between two patches in terms of parametric information and control points.
    bool IsEqual(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsCompatible(rOtherPatch))
            return false;

        // TODO compare the control points

        return true;
    }

    /// Compare between two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsEqual(rOtherPatch))
            return false;

        // TODO compare the grid function values

        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<TDim>& rOther)
    {
        return this->IsEqual(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Id = " << Id() << ", Add = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Neighbors = ";
        if (TDim == 2)
        {
            if (pLeft() != NULL)
                rOStream << " left:" << pLeft()->Id();
            if (pRight() != NULL)
                rOStream << " right:" << pRight()->Id();
            if (pTop() != NULL)
                rOStream << " top:" << pTop()->Id();
            if (pBottom() != NULL)
                rOStream << " bottom:" << pBottom()->Id();
        }
        else if (TDim == 3)
        {
            if (pLeft() != NULL)
                rOStream << " left:" << pLeft()->Id();
            if (pRight() != NULL)
                rOStream << " right:" << pRight()->Id();
            if (pTop() != NULL)
                rOStream << " top:" << pTop()->Id();
            if (pBottom() != NULL)
                rOStream << " bottom:" << pBottom()->Id();
            if (pFront() != NULL)
                rOStream << " front:" << pFront()->Id();
            if (pBack() != NULL)
                rOStream << " back:" << pBack()->Id();
        }
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
    typename Patch<TDim>::WeakPointer mpLeft;
    typename Patch<TDim>::WeakPointer mpRight;
    typename Patch<TDim>::WeakPointer mpTop;
    typename Patch<TDim>::WeakPointer mpBottom;
    typename Patch<TDim>::WeakPointer mpFront;
    typename Patch<TDim>::WeakPointer mpBack;

    /**
     * pointer to parent multipatch
     */
    typename MultiPatch<TDim>::WeakPointer mpParentMultiPatch;

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    /// Auxiliary
    template<class TGridFunctionType>
    void CheckSize(const TGridFunctionType& rGrid) const
    {
        if (rGrid.Size() != this->TotalNumber())
        {
            KRATOS_WATCH(rGrid.Size())
            KRATOS_WATCH(this->TotalNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The size of grid function is not compatible", "")
        }
    }
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
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
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
    rOStream << "-------------Begin PatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End PatchInfo-------------" << std::endl;
    return rOStream;
}


/**
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<std::size_t TDim>
class MultiPatch : public boost::enable_shared_from_this<MultiPatch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);

    /// Type definition
    typedef Patch<TDim> PatchType;
    typedef PointerVectorSet<PatchType, IndexedObject> PatchContainerType;

    /// Default constructor
    MultiPatch() : mGridSystemSize(0) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPatch(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
    }

    /// Reset Id for all the patches
    void ResetId()
    {
        std::size_t Id = 0;
        for (typename PatchContainerType::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->SetId(++Id);
        }
    }

    /// Validate the MultiPatch
    virtual bool Validate() const
    {
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            bool check = it->Validate();
            if (!check)
                return false;
        }

        return true;
    }

    /// Get the grid system size
    const std::size_t& GridSystemSize() const {return mGridSystemSize;}

    /// iterators
    typename PatchContainerType::iterator begin() {return mpPatches.begin();}
    typename PatchContainerType::const_iterator begin() const {return mpPatches.begin();}
    typename PatchContainerType::iterator end() {return mpPatches.end();}
    typename PatchContainerType::const_iterator end() const {return mpPatches.end();}

    /// Get the patch with specific Id
    typename PatchType::Pointer GetPatch(const std::size_t& Id) {return mpPatches(Id);}

    /// Access the underlying list of patches
    PatchContainerType& Patches() {return mpPatches;}

    /// Access the underlying list of patches
    const PatchContainerType& Patches() const {return mpPatches;}

    /// Get the number of patches
    const std::size_t& size() const {return mpPatches.size();}

    /// Enumerate all the patches
    void Enumerate(const std::size_t& start)
    {
        // firstly assign all the index to all the patch/grid functions to -1
        for (typename PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            (*it)->ResetDofs();
        }

        // secondly enumerate each patch
        mGridSystemSize = 0;
        std::set<std::size_t> enumerated_patches;
        for (typename PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            // enumerate the patch and remember
            mGridSystemSize = (*it)->Enumerate(mGridSystemSize);
            enumerated_patches.insert((*it)->Id());

            // for each patch, enumerate its neighbors
            for (typename PatchType::NeighborPatchContainerType::iterator it2 = (*it)->GetNeighbors().begin();
                    it2 != (*it)->GetNeighbors().end(); ++it2)
            {
                if ( std::find(enumerated_patches.begin(), enumerated_patches.end(), (*it2)->Id()) == enumerated_patches.end() )
                {
                    (*it2)->Enumerate(mGridSystemSize);
                    enumerated_patches.insert((*it2)->Id());
                }
            }
        }
    }

    /// Make the two patches neighbor. This required that two patches are conformed at the interface.
    /// If two patches are conformed, then the grid function on side1 will be transferred to side2.
    static void MakeNeighbor(typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2)
    {
        typename Patch<TDim-1>::Pointer pBPatch1 = pPatch1->ConstructBoundaryPatch(side1);
        typename Patch<TDim-1>::Pointer pBPatch2 = pPatch2->ConstructBoundaryPatch(side2);

        if( (*pBPatch1) == (*pBPatch2) )
        {
            // synchronize grid function data
            // temporarily disable. I think it can create potential conflict.
            // pBPatch1->SynchronizeGridFunction(*pBPatch2);

            // set the neighbor information
            if (side1 == _LEFT_)
                pPatch1->pSetLeft(pPatch2);
            else if (side1 == _RIGHT_)
                pPatch1->pSetRight(pPatch2);
            else if (side1 == _TOP_)
                pPatch1->pSetTop(pPatch2);
            else if (side1 == _BOTTOM_)
                pPatch1->pSetBottom(pPatch2);
            else if (side1 == _FRONT_)
                pPatch1->pSetFront(pPatch2);
            else if (side1 == _BACK_)
                pPatch1->pSetBack(pPatch2);

            if (side2 == _LEFT_)
                pPatch2->pSetLeft(pPatch1);
            else if (side2 == _RIGHT_)
                pPatch2->pSetRight(pPatch1);
            else if (side2 == _TOP_)
                pPatch2->pSetTop(pPatch1);
            else if (side2 == _BOTTOM_)
                pPatch2->pSetBottom(pPatch1);
            else if (side2 == _FRONT_)
                pPatch2->pSetFront(pPatch1);
            else if (side2 == _BACK_)
                pPatch2->pSetBack(pPatch1);

            // KRATOS_WATCH(*pPatch1)
            // KRATOS_WATCH(*pPatch2)
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The two patch's boundaries are not conformed", "")
    }

    /// Information

    void PrintAddress(typename PatchType::Pointer pPatch)
    {
        std::cout << pPatch << std::endl;
    }

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch overview: Number of patches = " << mpPatches.size();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch details:" << std::endl;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
            rOStream << (*it) << std::endl;
    }

protected:

    // /// Access the underlying list of patches
    // PatchContainerType& Patches() {return mpPatches;}

private:

    std::size_t mGridSystemSize;
    PatchContainerType mpPatches; // container for all the patches

};

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatch<TDim>& rThis)
{
    rOStream << ">>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    rOStream << "-------------Begin MultiPatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << "-------------End MultiPatchInfo-------------" << std::endl;
    rOStream << ">>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    return rOStream;
}



} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED defined


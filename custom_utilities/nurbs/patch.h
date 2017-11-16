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
#include "custom_utilities/nurbs/grid_function.h"


#define DEBUG_DESTROY

namespace Kratos
{


// Declaration
template<int TDim> class Patch;
template<int TDim> class MultiPatch;

/// Define the location of the sub-patches on the patch. The location for corners, edges or faces can be seen in
/// Fig.2, Burstedde et al, SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
enum SubPatchLocator
{
    _C0_, _C1_, _C2_, _C3_, _C4_, _C5_, _C6_, _C7_,
    _E0_, _E1_, _E2_, _E3_, _E4_, _E5_, _E6_, _E7_, _E8_, _E9_, _E10_, _E11_,
    _F0_, _F1_, _F2_, _F3_, _F4_, _F5_
};

/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<int TDim>
class Patch : public boost::enable_shared_from_this<Patch<TDim> >
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

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id), mFESpace(NULL)
    {
        mpNeighbors.resize(2*TDim);
        mpSubPatches.resize(TDim*(TDim-1)*(TDim-1)/2);
        mpSubSubPatches.resize(2*TDim*(TDim-1));
        mpSubSubSubPatches.resize(2^TDim);
    }

    /// Constructor with id and FESpace
    Patch(const std::size_t& Id, typename FESpace<TDim>::Pointer FESpace) : mId(Id), mFESpace(FESpace)
    {
        mpNeighbors.resize(2*TDim);
        mpSubPatches.resize(TDim*(TDim-1)*(TDim-1)/2);
        mpSubSubPatches.resize(2*TDim*(TDim-1));
        mpSubSubSubPatches.resize(2^TDim);
        if (mFESpace == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Invalid FESpace is provided", "")
    }

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

    /// Set the corresponding FESpace for the patch
    void SetFESpace(typename FESpace<TDim>::Pointer FESpace) {mFESpace = FESpace;}

    /// Get the FESpace
    typename FESpace<TDim>::Pointer FESpace() {return mFESpace;}
    typename FESpace<TDim>::ConstPointer FESpace() const {return mFESpace;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        assert(mFESpace == NULL);
        return mFESpace->TotalNumber();
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        assert(mFESpace == NULL);
        return mFESpace->Order(i);
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
    void ResetIds()
    {
        mFESpace->ResetIds();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Set the control point grid
    typename GridFunction<TDim, ControlPointType>::Pointer CreateControlPointGridFunction(typename ControlGrid<ControlPointType>::Pointer pControlPointGrid)
    {
        CheckSize(*pControlPointGrid, __FUNCTION__);
        mpControlPointGridFunc = typename GridFunction<TDim, ControlPointType>::Pointer(new GridFunction<TDim, ControlPointType>(mFESpace, pControlPointGrid));
        return mpControlPointGridFunc;
    }

    /// Get the control point grid function
    typename GridFunction<TDim, ControlPointType>::Pointer ControlPointGridFunction() {return mpControlPointGridFunc;}

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::ConstPointer ControlPointGridFunction() const {return mpControlPointGridFunc;}

    /// Get the control point weights vector
    std::vector<double> GetControlWeights() const
    {
        const typename GridFunction<TDim, ControlPointType>::DataContainerType& GridData = ControlPointGridFunction()->ControlGrid()->Data();
        std::vector<double> Weights(GridData.size());
        for (std::size_t i = 0; i < GridData.size(); ++i)
            Weights[i] = GridData[i].W();
        return Weights;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function for double variable
    typename DoubleGridFunctionType::Pointer CreateDoubleGridFunction(typename ControlGrid<double>::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename DoubleGridFunctionType::Pointer pNewGridFunc = typename DoubleGridFunctionType::Pointer(new DoubleGridFunctionType(mFESpace, pControlGrid));
        mpDoubleGridFuncs.push_back(pNewGridFunc);
        return pNewGridFunc;
    }

    /// Get the underlying double grid functions
    DoubleGridFunctionContainterType& DoubleGridFunctions() {return mpDoubleGridFuncs;}

    /// Get the underlying double grid functions
    const DoubleGridFunctionContainterType& DoubleGridFunctions() const {return mpDoubleGridFuncs;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function for array_1d variable
    typename Array1DGridFunctionType::Pointer CreateArray1DGridFunction(typename ControlGrid<array_1d<double, 3> >::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename Array1DGridFunctionType::Pointer pNewGridFunc = typename Array1DGridFunctionType::Pointer(new Array1DGridFunctionType(mFESpace, pControlGrid));
        mpArray1DGridFuncs.push_back(pNewGridFunc);
        return pNewGridFunc;
    }

    /// Get the underlying array_1d grid functions
    Array1DGridFunctionContainerType& Array1DGridFunctions() {return mpArray1DGridFuncs;}

    /// Get the underlying array_1d grid functions
    const Array1DGridFunctionContainerType& Array1DGridFunctions() const {return mpArray1DGridFuncs;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function for vector variable
    typename VectorGridFunctionType::Pointer CreateVectorGridFunction(typename ControlGrid<Vector>::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename VectorGridFunctionType::Pointer pNewGridFunc = typename VectorGridFunctionType::Pointer(new VectorGridFunctionType(mFESpace, pControlGrid));
        mpVectorGridFuncs.push_back(pNewGridFunc);
        return pNewGridFunc;
    }

    /// Get the underlying vector grid functions
    VectorGridFunctionContainerType& VectorGridFunctions() {return mpVectorGridFuncs;}

    /// Get the underlying vector grid functions
    const VectorGridFunctionContainerType& VectorGridFunctions() const {return mpVectorGridFuncs;}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Validate the patch
    virtual bool Validate() const
    {
        if (Id() == 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The patch must have an Id", "")
        }

        if (mpControlPointGridFunc->ControlGrid()->Size() != this->TotalNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        for (typename DoubleGridFunctionContainterType::const_iterator it = DoubleGridFunctions().begin();
                it != DoubleGridFunctions().end(); ++it)
        {
            if ((*it)->ControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The double variable grid is incompatible", (*it)->ControlGrid()->Name())
                return false;
            }
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions().begin();
                it != Array1DGridFunctions().end(); ++it)
        {
            if ((*it)->ControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The array_1d variable grid is incompatible", (*it)->ControlGrid()->Name())
                return false;
            }
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions().begin();
                it != VectorGridFunctions().end(); ++it)
        {
            if ((*it)->ControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The vector variable grid is incompatible", (*it)->ControlGrid()->Name())
                return false;
            }
        }

        // check the compatibility between patch
        if (pNeighbor(_LEFT_) != NULL)
        {
            if (this->Type() != pNeighbor(_LEFT_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and left neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _LEFT_, *pNeighbor(_LEFT_), _RIGHT_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and left neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pNeighbor(_RIGHT_) != NULL)
        {
            if (this->Type() != pNeighbor(_RIGHT_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and right neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _RIGHT_, *pNeighbor(_LEFT_), _LEFT_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and right neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pNeighbor(_TOP_) != NULL)
        {
            if (this->Type() != pNeighbor(_TOP_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and right neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _TOP_, *pNeighbor(_TOP_), _BOTTOM_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and top neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pNeighbor(_BOTTOM_) != NULL)
        {
            if (this->Type() != pNeighbor(_BOTTOM_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and bottom neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _BOTTOM_, *pNeighbor(_TOP_), _TOP_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and bottom neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pNeighbor(_FRONT_) != NULL)
        {
            if (this->Type() != pNeighbor(_FRONT_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and front neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _FRONT_, *pNeighbor(_FRONT_), _BACK_);
                if (!check)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and front neighbor is incompatible", "")
                    return false;
                }
            }
        }

        if (pNeighbor(_BACK_) != NULL)
        {
            if (this->Type() != pNeighbor(_BACK_)->Type())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and back neighbor is incompatible", "")
            }
            else
            {
                bool check = CheckBoundaryCompatibility(*this, _BACK_, *pNeighbor(_BACK_), _FRONT_);
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
        typename Patch<TDim-1>::Pointer pBPatch = typename Patch<TDim-1>::Pointer(new Patch<TDim-1>(-1));

        typename FESpace<TDim-1>::Pointer pBFESpace = this->FESpace()->ConstructBoundaryFESpace(side);
        pBPatch->SetFESpace(pBFESpace);

        // TODO transfer the control values

        return pBPatch;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    void pSetNeighbor(const BoundarySide& side, typename Patch<TDim>::Pointer pNeighbor) {mpNeighbors[side] = pNeighbor->shared_from_this();}
    Patch<TDim>& Neighbor(const BoundarySide& side) {return *pNeighbor(side);}
    const Patch<TDim>& Neighbor(const BoundarySide& side) const {return *pNeighbor(side);}
    typename Patch<TDim>::Pointer pNeighbor(const BoundarySide& side)
    {
        if (side < 2*TDim)
            return mpNeighbors[side].lock();
        else
            return NULL;
    }
    typename Patch<TDim>::ConstPointer pNeighbor(const BoundarySide& side) const
    {
        if (side < 2*TDim)
            return mpNeighbors[side].lock();
        else
            return NULL;
    }

    /// Get/Set the parent multipatch
    void pSetParentMultiPatch(typename MultiPatch<TDim>::Pointer pPatch) {mpParentMultiPatch = pPatch;}
    MultiPatch<TDim>& ParentMultiPatch() {return *pParentMultiPatch();}
    const MultiPatch<TDim>& ParentMultiPatch() const {return *pParentMultiPatch();}
    typename MultiPatch<TDim>::Pointer pParentMultiPatch() {return mpParentMultiPatch.lock();}
    const typename MultiPatch<TDim>::Pointer pParentMultiPatch() const {return mpParentMultiPatch.lock();}

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two patches in terms of its parametric information. The grid function data, including control points, shall not be checked.
    virtual bool IsCompatible(const Patch<TDim>& rOtherPatch) const
    {
        return *(this->FESpace()) == *(rOtherPatch.FESpace());
    }

    /// Compare between two patches in terms of parametric information and control points.
    bool IsEquivalent(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsCompatible(rOtherPatch))
            return false;

        // TODO compare the control points

        return true;
    }

    /// Compare between two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsEquivalent(rOtherPatch))
            return false;

        // TODO compare the grid function values

        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<TDim>& rOther)
    {
        return (Id() == rOther.Id()) && this->IsSame(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Id = " << Id() << ", Add = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        if (FESpace() != NULL)
            rOStream << *FESpace() << std::endl;
        if (ControlPointGridFunction() != NULL)
            rOStream << *(ControlPointGridFunction()->ControlGrid());
        rOStream << "Neighbors = ";
        if (TDim == 2)
        {
            if (pNeighbor(_LEFT_) != NULL)
                rOStream << " left:" << pNeighbor(_LEFT_)->Id();
            if (pNeighbor(_RIGHT_) != NULL)
                rOStream << " right:" << pNeighbor(_RIGHT_)->Id();
            if (pNeighbor(_TOP_) != NULL)
                rOStream << " top:" << pNeighbor(_TOP_)->Id();
            if (pNeighbor(_BOTTOM_) != NULL)
                rOStream << " bottom:" << pNeighbor(_BOTTOM_)->Id();
        }
        else if (TDim == 3)
        {
            if (pNeighbor(_LEFT_) != NULL)
                rOStream << " left:" << pNeighbor(_LEFT_)->Id();
            if (pNeighbor(_RIGHT_) != NULL)
                rOStream << " right:" << pNeighbor(_RIGHT_)->Id();
            if (pNeighbor(_TOP_) != NULL)
                rOStream << " top:" << pNeighbor(_TOP_)->Id();
            if (pNeighbor(_BOTTOM_) != NULL)
                rOStream << " bottom:" << pNeighbor(_BOTTOM_)->Id();
            if (pNeighbor(_FRONT_) != NULL)
                rOStream << " front:" << pNeighbor(_FRONT_)->Id();
            if (pNeighbor(_BACK_) != NULL)
                rOStream << " back:" << pNeighbor(_BACK_)->Id();
        }
    }

private:

    std::size_t mId;

    // shape function information
    typename FESpace<TDim>::Pointer mFESpace;

    // grid functions, including the grid function over the control points
    typename GridFunction<TDim, ControlPointType>::Pointer mpControlPointGridFunc;
    DoubleGridFunctionContainterType mpDoubleGridFuncs;
    Array1DGridFunctionContainerType mpArray1DGridFuncs;
    VectorGridFunctionContainerType mpVectorGridFuncs;

    /**
     * neighboring data
     */
    std::vector<typename Patch<TDim>::WeakPointer> mpNeighbors;

    /**
     * sub-patch data
     */
    std::vector<typename Patch<TDim-1>::WeakPointer> mpSubPatches;
    std::vector<typename Patch<TDim-2>::WeakPointer> mpSubSubPatches;
    std::vector<typename Patch<TDim-3>::WeakPointer> mpSubSubSubPatches;

    /**
     * pointer to parent multipatch
     */
    typename MultiPatch<TDim>::WeakPointer mpParentMultiPatch;

    /// Empty Constructor for serializer
    Patch() : mId(0), mFESpace(NULL) {}

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
    void CheckSize(const TGridFunctionType& rGrid, const std::string& source) const
    {
        if (rGrid.Size() != this->TotalNumber())
        {
            std::stringstream ss;
            ss << "The size of grid function (" << rGrid.size() << ") is not compatible with the current number of control values (" << this->TotalNumber()
               << ") of patch " << Id() << ". Error at " << source;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }
};



/**
 * Template specific instantiation for null-D patch to terminate the compilation.
 * In fact, null-D patch is a vertex
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

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<0>::Pointer FESpace) {}

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
        return Id() == rOther.Id();
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

/**
 * Template specific instantiation for -1-D patch to terminate the compilation.
 */
template<>
class Patch<-1>
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

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<-1>::Pointer FESpace) {}

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
        return "Patch<-1>D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<-1>& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    virtual bool CheckBoundaryCompatibility(const Patch<-1>& rPatch1, const BoundarySide& side1,
            const Patch<-1>& rPatch2, const BoundarySide& side2) const
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
        rOStream << "Patch<-1>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/**
 * Template specific instantiation for -2-D patch to terminate the compilation.
 */
template<>
class Patch<-2>
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

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<-2>::Pointer FESpace) {}

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
        return "Patch<-2>D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<-2>& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    virtual bool CheckBoundaryCompatibility(const Patch<-2>& rPatch1, const BoundarySide& side1,
            const Patch<-2>& rPatch2, const BoundarySide& side2) const
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-2>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<-2>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const Patch<TDim>& rThis)
{
    rOStream << "-------------Begin PatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End PatchInfo-------------";
    return rOStream;
}


/**
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<int TDim>
class MultiPatch : public boost::enable_shared_from_this<MultiPatch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);

    /// Type definition
    typedef Patch<TDim> PatchType;
    typedef PointerVectorSet<PatchType, IndexedObject> PatchContainerType;

    typedef Patch<0> VertexPatchType;
    typedef PointerVectorSet<VertexPatchType, IndexedObject> VertexPatchContainerType;

    typedef Patch<1> EdgePatchType;
    typedef PointerVectorSet<EdgePatchType, IndexedObject> EdgePatchContainerType;

    typedef Patch<2> FacePatchType;
    typedef PointerVectorSet<FacePatchType, IndexedObject> FacePatchContainerType;

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
    /// WARNING!!! be careful with this routine
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
        // if( pBPatch1->IsCompatible(*pBPatch2) )
        {
            // synchronize grid function data
            // temporarily disable. I think it can create potential conflict.
            // pBPatch1->SynchronizeGridFunction(*pBPatch2);

            // set the neighbor information
            pPatch1->pSetNeighbor(side1, pPatch2);
            pPatch2->pSetNeighbor(side2, pPatch1);

            // KRATOS_WATCH(*pPatch1)
            // KRATOS_WATCH(*pPatch2)

            // TODO create the multipatch topology structure
            
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The two patch's boundaries are not conformed", "")
    }

    /// Generate the multipatch topology that can be read in Glvis
    void GenerateCornerTopology(std::size_t& nvertices,
        std::vector<std::vector<std::size_t> >& elements,
        std::vector<std::vector<std::size_t> >& boundary,
        std::vector<std::size_t>& boundary_attr,
        std::vector<std::size_t>& edges_knotv) const
    {

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

private:

    std::size_t mGridSystemSize;
    PatchContainerType mpPatches; // container for all the patches

    VertexPatchContainerType mpVertexPatches;
    EdgePatchContainerType mpEdgePatches;
    FacePatchContainerType mpFacePatches;
};

/// output stream function
template<int TDim>
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


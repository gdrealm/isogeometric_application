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
#include <tuple>

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

    typedef std::size_t vertex_t;
    typedef std::tuple<std::size_t, std::size_t, std::size_t, int> edge_t;
    //                  vertex1     vertex2     knot index   is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, int> face_t;
    //                  vertex1     vertex2     vertex3         vertex4     is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t> volume_t;


    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id), mFESpace(NULL)
    {
        mpNeighbors.resize(2*TDim);
    }

    /// Constructor with id and FESpace
    Patch(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace) : mId(Id), mFESpace(pFESpace)
    {
        mpNeighbors.resize(2*TDim);
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
    void SetFESpace(typename FESpace<TDim>::Pointer pFESpace) {mFESpace = pFESpace;}

    /// Get the FESpace pointer
    typename FESpace<TDim>::Pointer pFESpace() {return mFESpace;}

    /// Get the FESpace pointer
    typename FESpace<TDim>::ConstPointer pFESpace() const {return mFESpace;}

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
        const typename GridFunction<TDim, ControlPointType>::DataContainerType& GridData = ControlPointGridFunction()->pControlGrid()->Data();
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

        if (mpControlPointGridFunc->pControlGrid()->Size() != this->TotalNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        for (typename DoubleGridFunctionContainterType::const_iterator it = DoubleGridFunctions().begin();
                it != DoubleGridFunctions().end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The double variable grid is incompatible", (*it)->pControlGrid()->Name())
                return false;
            }
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions().begin();
                it != Array1DGridFunctions().end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The array_1d variable grid is incompatible", (*it)->pControlGrid()->Name())
                return false;
            }
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions().begin();
                it != VectorGridFunctions().end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The vector variable grid is incompatible", (*it)->pControlGrid()->Name())
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
    static bool CheckBoundaryCompatibility(const Patch<TDim>& rPatch1, const BoundarySide& side1,
            const Patch<TDim>& rPatch2, const BoundarySide& side2)
    {
        typename Patch<TDim-1>::Pointer BPatch1 = rPatch1.ConstructBoundaryPatch(side1);
        typename Patch<TDim-1>::Pointer BPatch2 = rPatch1.ConstructBoundaryPatch(side2);

        return (*BPatch1) == (*BPatch2);
    }

    /// Check the boundary compatibility between this patch and the other patch
    bool CheckBoundaryCompatibility(const BoundarySide& side1,
            const Patch<TDim>& rOtherPatch, const BoundarySide& side2) const
    {
        return CheckBoundaryCompatibility(*this, side1, rOtherPatch, side2);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Construct the boundary patch based on side
    virtual typename Patch<TDim-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    {
        typename Patch<TDim-1>::Pointer pBPatch = typename Patch<TDim-1>::Pointer(new Patch<TDim-1>(-1));

        typename FESpace<TDim-1>::Pointer pBFESpace = this->pFESpace()->ConstructBoundaryFESpace(side);
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

    /// Generate topology data to visualize with GLVis
    void GenerateTopolgyData(std::size_t& starting_vertex_id,
            std::vector<vertex_t>& vertices,
            std::vector<edge_t>& edges,
            std::vector<face_t>& faces,
            std::vector<volume_t>& volumes,
            std::size_t& starting_knotv_id,
            std::vector<std::size_t>& knotv ) const
    {
        if (TDim == 1)
        {
            vertices.resize(2);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;

            knotv.resize(0);
            knotv[0] = starting_knotv_id++;

            edges.resize(1);
            edges[0] = std::make_tuple(vertices[0], vertices[1], knotv[0], 0);

            faces.resize(0);
            volumes.resize(0);
        }
        else if (TDim == 2)
        {
            /// Reference for edge mapping: Fig.2, Burstedde et al, SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
            /// Mapping for edges: table 2
            vertices.resize(4);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;
            vertices[2] = starting_vertex_id++;
            vertices[3] = starting_vertex_id++;

            knotv.resize(2);
            knotv[0] = starting_knotv_id++;
            knotv[1] = starting_knotv_id++;

            edges.resize(4);
            edges[0] = std::make_tuple(vertices[0], vertices[2], knotv[1], 1);
            edges[1] = std::make_tuple(vertices[1], vertices[3], knotv[1], 1);
            edges[2] = std::make_tuple(vertices[0], vertices[1], knotv[0], 1);
            edges[3] = std::make_tuple(vertices[2], vertices[3], knotv[0], 1);

            faces.resize(1);
            faces[0] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], 0);

            volumes.resize(0);
        }
        else if (TDim == 3)
        {
            /// Reference for edge mapping: Fig.2, Burstedde et al, SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
            /// Mapping for faces: table 2
            vertices.resize(8);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;
            vertices[2] = starting_vertex_id++;
            vertices[3] = starting_vertex_id++;
            vertices[4] = starting_vertex_id++;
            vertices[5] = starting_vertex_id++;
            vertices[6] = starting_vertex_id++;
            vertices[7] = starting_vertex_id++;

            knotv.resize(3);
            knotv[0] = starting_knotv_id++;
            knotv[1] = starting_knotv_id++;
            knotv[2] = starting_knotv_id++;

            edges.resize(12);
            edges[0] = std::make_tuple(vertices[0], vertices[1], knotv[0], 1);
            edges[1] = std::make_tuple(vertices[2], vertices[3], knotv[0], 1);
            edges[2] = std::make_tuple(vertices[4], vertices[5], knotv[0], 1);
            edges[3] = std::make_tuple(vertices[6], vertices[7], knotv[0], 1);
            edges[4] = std::make_tuple(vertices[0], vertices[2], knotv[1], 1);
            edges[5] = std::make_tuple(vertices[1], vertices[3], knotv[1], 1);
            edges[6] = std::make_tuple(vertices[4], vertices[6], knotv[1], 1);
            edges[7] = std::make_tuple(vertices[5], vertices[7], knotv[1], 1);
            edges[8] = std::make_tuple(vertices[0], vertices[4], knotv[2], 1);
            edges[9] = std::make_tuple(vertices[1], vertices[5], knotv[2], 1);
            edges[10] = std::make_tuple(vertices[2], vertices[6], knotv[2], 1);
            edges[11] = std::make_tuple(vertices[3], vertices[7], knotv[2], 1);

            faces.resize(6);
            faces[0] = std::make_tuple(vertices[0], vertices[2], vertices[4], vertices[6], 1);
            faces[1] = std::make_tuple(vertices[1], vertices[3], vertices[5], vertices[7], 1);
            faces[2] = std::make_tuple(vertices[0], vertices[1], vertices[4], vertices[5], 1);
            faces[3] = std::make_tuple(vertices[2], vertices[3], vertices[6], vertices[7], 1);
            faces[4] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], 1);
            faces[5] = std::make_tuple(vertices[4], vertices[5], vertices[6], vertices[7], 1);

            volumes.resize(1);
            volumes[0] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5], vertices[6], vertices[7]);
        }
        else
        {
            std::stringstream ss;
            ss << __FUNCTION__ << " is not implemented for " << TDim;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two patches in terms of its parametric information. The grid function data, including control points, shall not be checked.
    virtual bool IsCompatible(const Patch<TDim>& rOtherPatch) const
    {
        return *(this->pFESpace()) == *(rOtherPatch.pFESpace());
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
        if (pFESpace() != NULL)
            rOStream << *pFESpace() << std::endl;
        if (ControlPointGridFunction() != NULL)
            rOStream << *(ControlPointGridFunction()->pControlGrid());
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
    void SetFESpace(typename FESpace<0>::Pointer pFESpace) {}

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
    static bool CheckBoundaryCompatibility(const Patch<0>& rPatch1, const BoundarySide& side1,
            const Patch<0>& rPatch2, const BoundarySide& side2)
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
    void SetFESpace(typename FESpace<-1>::Pointer pFESpace) {}

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
    static bool CheckBoundaryCompatibility(const Patch<-1>& rPatch1, const BoundarySide& side1,
            const Patch<-1>& rPatch2, const BoundarySide& side2)
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
    void SetFESpace(typename FESpace<-2>::Pointer pFESpace) {}

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
    static bool CheckBoundaryCompatibility(const Patch<-2>& rPatch1, const BoundarySide& side1,
            const Patch<-2>& rPatch2, const BoundarySide& side2)
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

    typedef typename Patch<TDim>::vertex_t vertex_t;
    typedef typename Patch<TDim>::edge_t edge_t;
    typedef typename Patch<TDim>::face_t face_t;
    typedef typename Patch<TDim>::volume_t volume_t;

    /// Default constructor
    MultiPatch() {}

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
    void Enumerate(std::size_t& EquationSystemSize, std::set<std::size_t>& enumerated_patches)
    {
        // secondly enumerate each patch
        EquationSystemSize = 0;
        for (typename PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            // enumerate the patch and remember
            if ( std::find(enumerated_patches.begin(), enumerated_patches.end(), it->Id()) == enumerated_patches.end() )
            {
                EquationSystemSize = it->pFESpace()->Enumerate(EquationSystemSize);
                enumerated_patches.insert(it->Id());
                // KRATOS_WATCH(EquationSystemSize)

                // find the neighbors and enumerate the boundary of the neighbors
                if (it->pNeighbor(_LEFT_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_LEFT_, *(it->pNeighbor(_LEFT_)), _RIGHT_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with LEFT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_LEFT_);
                        it->pNeighbor(_LEFT_)->pFESpace()->AssignBoundaryFunctionIndices(_RIGHT_, func_indices);
                    }
                }

                if (it->pNeighbor(_RIGHT_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_RIGHT_, *(it->pNeighbor(_RIGHT_)), _LEFT_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with RIGHT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_RIGHT_);
                        it->pNeighbor(_RIGHT_)->pFESpace()->AssignBoundaryFunctionIndices(_LEFT_, func_indices);
                    }
                }

                if (it->pNeighbor(_TOP_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_TOP_, *(it->pNeighbor(_TOP_)), _BOTTOM_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with RIGHT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_TOP_);
                        it->pNeighbor(_TOP_)->pFESpace()->AssignBoundaryFunctionIndices(_BOTTOM_, func_indices);
                    }
                }

                if (it->pNeighbor(_BOTTOM_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_BOTTOM_, *(it->pNeighbor(_BOTTOM_)), _TOP_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with RIGHT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_BOTTOM_);
                        it->pNeighbor(_BOTTOM_)->pFESpace()->AssignBoundaryFunctionIndices(_TOP_, func_indices);
                    }
                }

                if (it->pNeighbor(_FRONT_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_FRONT_, *(it->pNeighbor(_FRONT_)), _BACK_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with RIGHT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_FRONT_);
                        it->pNeighbor(_FRONT_)->pFESpace()->AssignBoundaryFunctionIndices(_BACK_, func_indices);
                    }
                }

                if (it->pNeighbor(_BACK_) != NULL)
                {
                    // check the boundary compatibility again
                    if (!it->CheckBoundaryCompatibility(_BACK_, *(it->pNeighbor(_BACK_)), _FRONT_))
                    {
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with RIGHT neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = it->pFESpace()->ExtractBoundaryFunctionIndices(_BACK_);
                        it->pNeighbor(_BACK_)->pFESpace()->AssignBoundaryFunctionIndices(_FRONT_, func_indices);
                    }
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
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The two patch's boundaries are not conformed", "")
    }

    /// Generate the multipatch topology that can be read in GLVis
    void GenerateCornerTopology(std::size_t& nvertices,
        std::vector<std::vector<std::size_t> >& elements,
        std::vector<std::vector<std::size_t> >& boundary,
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t, int> >& edges,
        std::map<std::size_t, std::vector<double> >& knotvecs) const
    {
        std::map<std::size_t, std::vector<vertex_t> > patch_vertices;
        std::map<std::size_t, std::vector<edge_t> > patch_edges;
        std::map<std::size_t, std::vector<face_t> > patch_faces;
        std::map<std::size_t, std::vector<volume_t> > patch_volumes;
        std::map<std::size_t, std::vector<std::size_t> > patch_knotv;

        // generate vertices, edges, faces for all the patches
        std::size_t start_vertex_id = 0;
        std::size_t start_knotv_id = 0;
        std::map<std::size_t, std::vector<double> > all_knotvec;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            it->GenerateTopolgyData(start_vertex_id, patch_vertices[id], patch_edges[id], patch_faces[id], patch_volumes[id], start_knotv_id, patch_knotv[id]);

            typename NURBSFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<NURBSFESpace<TDim> >(it->pFESpace());

            // collect the knot vector in respective dimension
            for (std::size_t i = 0; i < TDim; ++i)
            {
                for (std::size_t j = 0; j < pFESpace->KnotVector(i).size(); ++j)
                    all_knotvec[patch_knotv[id][i]].push_back(pFESpace->KnotVector(i)[j]);
            }
        }

        #ifdef DEBUG_GLVIS_EXPORT
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();

            std::cout << "edge (p1) for patch " << id << std::endl;
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::cout << std::get<2>(patch_edges[id][i])
                          << " " << std::get<0>(patch_edges[id][i])
                          << " " << std::get<1>(patch_edges[id][i])
                          << " " << std::get<3>(patch_edges[id][i])
                          << std::endl;
            }
        }
        #endif

        // for all patch, account for the corners and then renumbering the vertex and edge
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();

            if (it->pNeighbor(_LEFT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_LEFT_)->Id();
                SynchronizeVertices(_LEFT_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][0]) = patch_knotv[id][1]; // change the knot index vector
                    std::get<2>(patch_edges[other_id][1]) = patch_knotv[id][1];
                    std::get<3>(patch_edges[other_id][1]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][0]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    patch_knotv[other_id][2] = patch_knotv[id][2];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_RIGHT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_RIGHT_)->Id();
                SynchronizeVertices(_RIGHT_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][0]) = patch_knotv[id][1];
                    std::get<2>(patch_edges[other_id][1]) = patch_knotv[id][1];
                    std::get<3>(patch_edges[other_id][0]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][1]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    patch_knotv[other_id][2] = patch_knotv[id][2];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_TOP_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_TOP_)->Id();
                SynchronizeVertices(_TOP_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][2]) = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][3]) = patch_knotv[id][0];
                    std::get<3>(patch_edges[other_id][2]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][3]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_BOTTOM_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_BOTTOM_)->Id();
                SynchronizeVertices(_BOTTOM_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                if (TDim == 2)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][2]) = patch_knotv[id][0];
                    std::get<2>(patch_edges[other_id][3]) = patch_knotv[id][0];
                    std::get<3>(patch_edges[other_id][3]) = 0; // disable the boundary flag
                    std::get<3>(patch_edges[id][2]) = 0; // disable the boundary flag
                }
                else if (TDim == 3)
                {
                    patch_knotv[other_id][0] = patch_knotv[id][0];
                    patch_knotv[other_id][1] = patch_knotv[id][1];
                    // TODO change knot vector index for edge
                    // TODO disable the boundary flag for edges and faces
                }
            }

            if (it->pNeighbor(_FRONT_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_FRONT_)->Id();
                SynchronizeVertices(_FRONT_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                patch_knotv[other_id][0] = patch_knotv[id][0];
                patch_knotv[other_id][2] = patch_knotv[id][2];
                // TODO change knot vector index for edge
                // TODO disable the boundary flag for edges and faces
            }

            if (it->pNeighbor(_BACK_) != NULL)
            {
                const std::size_t& other_id = it->pNeighbor(_BACK_)->Id();
                SynchronizeVertices(_BACK_, patch_vertices[id], patch_vertices[other_id], patch_edges[other_id], patch_faces[other_id], patch_volumes[other_id]);
                patch_knotv[other_id][0] = patch_knotv[id][0];
                patch_knotv[other_id][2] = patch_knotv[id][2];
                // TODO change knot vector index for edge
                // TODO disable the boundary flag for edges and faces
            }
        }

        #ifdef DEBUG_GLVIS_EXPORT
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();

            std::cout << "edge (p2) for patch " << id << std::endl;
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::cout << std::get<2>(patch_edges[id][i])
                          << " " << std::get<0>(patch_edges[id][i])
                          << " " << std::get<1>(patch_edges[id][i])
                          << " " << std::get<3>(patch_edges[id][i])
                          << std::endl;
            }
        }
        #endif

        // collect all the vertices in the multipatch
        std::set<vertex_t> all_vertices;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            all_vertices.insert(patch_vertices[id].begin(), patch_vertices[id].end());
        }

        // reassign each id a new index
        std::map<vertex_t, vertex_t> old_to_new;
        start_vertex_id = 0;
        for (typename std::set<vertex_t>::iterator it = all_vertices.begin(); it != all_vertices.end(); ++it)
            old_to_new[*it] = start_vertex_id++;

        // finally reassign the new id for all patches
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();

            for (std::size_t i = 0; i < patch_vertices[id].size(); ++i)
            {
                patch_vertices[id][i] = old_to_new[patch_vertices[id][i]];
            }

            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::get<0>(patch_edges[id][i]) = old_to_new[std::get<0>(patch_edges[id][i])];
                std::get<1>(patch_edges[id][i]) = old_to_new[std::get<1>(patch_edges[id][i])];
            }

            for (std::size_t i = 0; i < patch_faces[id].size(); ++i)
            {
                std::get<0>(patch_faces[id][i]) = old_to_new[std::get<0>(patch_faces[id][i])];
                std::get<1>(patch_faces[id][i]) = old_to_new[std::get<1>(patch_faces[id][i])];
                std::get<2>(patch_faces[id][i]) = old_to_new[std::get<2>(patch_faces[id][i])];
                std::get<3>(patch_faces[id][i]) = old_to_new[std::get<3>(patch_faces[id][i])];
            }

            for (std::size_t i = 0; i < patch_volumes[id].size(); ++i)
            {
                std::get<0>(patch_volumes[id][i]) = old_to_new[std::get<0>(patch_volumes[id][i])];
                std::get<1>(patch_volumes[id][i]) = old_to_new[std::get<1>(patch_volumes[id][i])];
                std::get<2>(patch_volumes[id][i]) = old_to_new[std::get<2>(patch_volumes[id][i])];
                std::get<3>(patch_volumes[id][i]) = old_to_new[std::get<3>(patch_volumes[id][i])];
                std::get<4>(patch_volumes[id][i]) = old_to_new[std::get<4>(patch_volumes[id][i])];
                std::get<5>(patch_volumes[id][i]) = old_to_new[std::get<5>(patch_volumes[id][i])];
                std::get<6>(patch_volumes[id][i]) = old_to_new[std::get<6>(patch_volumes[id][i])];
                std::get<7>(patch_volumes[id][i]) = old_to_new[std::get<7>(patch_volumes[id][i])];
            }
        }

        // collect all the knot vector index in all the patches
        std::set<std::size_t> all_knotvs;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                all_knotvs.insert(patch_knotv[id][i]);
        }

        #ifdef DEBUG_GLVIS_EXPORT
        std::cout << "all_knotvs:";
        for (std::set<std::size_t>::iterator it = all_knotvs.begin(); it != all_knotvs.end(); ++it)
            std::cout << " " << *it;
        std::cout << std::endl;
        #endif

        // reassign each knot vector a new index
        std::map<std::size_t, std::size_t> old_to_new_knotv;
        start_knotv_id = 0;
        for (std::set<std::size_t>::iterator it = all_knotvs.begin(); it != all_knotvs.end(); ++it)
            old_to_new_knotv[*it] = start_knotv_id++;

        // finally reassign the new id for all patches
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                patch_knotv[id][i] = old_to_new_knotv[patch_knotv[id][i]];
        }

        // reassign the knot vector index in each edge
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                std::get<2>(patch_edges[id][i]) = old_to_new_knotv[std::get<2>(patch_edges[id][i])];
            }
        }

        // assign the knot vector accordingly
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < TDim; ++i)
                knotvecs[ old_to_new_knotv[patch_knotv[id][i]] ] = all_knotvec[ patch_knotv[id][i] ];
        }

        // collect all the edges in all the patches
        std::set<edge_t> all_edges;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            const std::size_t& id = it->Id();
            for (std::size_t i = 0; i < patch_edges[id].size(); ++i)
            {
                all_edges.insert(patch_edges[id][i]);
            }
        }
        edges.assign(all_edges.begin(), all_edges.end());

        ///////////////////////////////////////////////////

        // export the information
        nvertices = all_vertices.size();

        if (TDim == 2)
        {
            for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
            {
                const std::size_t& id = it->Id();

                std::vector<std::size_t> elem(4);
                elem[0] = std::get<0>(patch_faces[id][0]);
                elem[1] = std::get<1>(patch_faces[id][0]);
                elem[2] = std::get<3>(patch_faces[id][0]); // here we switch the role of the vertex because GLVis only accepts the sequence 0-1-3-2 for quadrilateral
                elem[3] = std::get<2>(patch_faces[id][0]);

                elements.push_back(elem);
            }
        }
        else if (TDim == 3)
        {

        }
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

    PatchContainerType mpPatches; // container for all the patches

    /// Synchronize from 1->2
    void SynchronizeVertices( const BoundarySide& side1,
        std::vector<vertex_t>& vertices1,
        std::vector<vertex_t>& vertices2,
        std::vector<edge_t>& edges2,
        std::vector<face_t>& faces2,
        std::vector<volume_t>& volumes2 ) const
    {
        if (vertices1.size() != vertices2.size())
            KRATOS_THROW_ERROR(std::logic_error, "The number of vertices is not compatible", "")

        std::vector<int> map = GetJointMapping(TDim, side1);

        std::map<std::size_t, std::size_t> old_to_new;
        for (std::size_t i = 0; i < vertices2.size(); ++i)
            old_to_new[vertices2[i]] = vertices2[i];
        for (std::size_t i = 0; i < map.size()/2; ++i)
        {
            old_to_new[vertices2[map[i*2+1]]] = vertices1[map[i*2]];
            vertices2[map[i*2+1]] = vertices1[map[i*2]];
        }

        for (std::size_t i = 0; i < edges2.size(); ++i)
        {
            std::get<0>(edges2[i]) = old_to_new[std::get<0>(edges2[i])];
            std::get<1>(edges2[i]) = old_to_new[std::get<1>(edges2[i])];
        }

        for (std::size_t i = 0; i < faces2.size(); ++i)
        {
            std::get<0>(faces2[i]) = old_to_new[std::get<0>(faces2[i])];
            std::get<1>(faces2[i]) = old_to_new[std::get<1>(faces2[i])];
            std::get<2>(faces2[i]) = old_to_new[std::get<2>(faces2[i])];
            std::get<3>(faces2[i]) = old_to_new[std::get<3>(faces2[i])];
        }

        for (std::size_t i = 0; i < volumes2.size(); ++i)
        {
            std::get<0>(volumes2[i]) = old_to_new[std::get<0>(volumes2[i])];
            std::get<1>(volumes2[i]) = old_to_new[std::get<1>(volumes2[i])];
            std::get<2>(volumes2[i]) = old_to_new[std::get<2>(volumes2[i])];
            std::get<3>(volumes2[i]) = old_to_new[std::get<3>(volumes2[i])];
            std::get<4>(volumes2[i]) = old_to_new[std::get<4>(volumes2[i])];
            std::get<5>(volumes2[i]) = old_to_new[std::get<5>(volumes2[i])];
            std::get<6>(volumes2[i]) = old_to_new[std::get<6>(volumes2[i])];
            std::get<7>(volumes2[i]) = old_to_new[std::get<7>(volumes2[i])];
        }
    }

    std::vector<int> GetJointMapping(const int& dim, const BoundarySide& side) const
    {
        if (dim == 2)
        {
            if (side == _LEFT_)        return std::vector<int>{ 0, 1, /**/ 2, 3 };
            else if (side == _RIGHT_)  return std::vector<int>{ 1, 0, /**/ 3, 2 };
            else if (side == _TOP_)    return std::vector<int>{ 2, 0, /**/ 3, 1 };
            else if (side == _BOTTOM_) return std::vector<int>{ 0, 2, /**/ 1, 3 };
        }
        else if (dim == 3)
        {
            // TODO
        }
    }
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


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
#include <boost/any.hpp>
#include <boost/array.hpp>
#include <boost/enable_shared_from_this.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/array_1d.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/control_grid_utility.h"
#include "isogeometric_application/isogeometric_application.h"


#define DEBUG_DESTROY

namespace Kratos
{


// Declaration
template<int TDim> class Patch;
template<int TDim> class MultiPatch;


enum IsogeometricEchoFlags
{
    ECHO_REFIMENT   = 0b0000000000000001,
};


/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
 */
template<int TDim>
class Patch : public boost::enable_shared_from_this<Patch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef Transformation<double> TransformationType;

    typedef GridFunction<TDim, double> DoubleGridFunctionType;
    typedef std::vector<typename DoubleGridFunctionType::Pointer> DoubleGridFunctionContainerType;

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

    typedef FESpace<TDim> FESpaceType;

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
        std::cout << Type() << ", Id = " << Id() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Helper function to create new patch pointer
    static typename Patch<TDim>::Pointer Create(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace));
    }

    /// Get the working space dimension of the patch
    std::size_t WorkingSpaceDimension() const {return TDim;}

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
        if (i >= TDim) return 0;
        else return mFESpace->Order(i);
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
        return this->CreateGridFunction(CONTROL_POINT, pControlPointGrid);
    }

    /// Get the control point grid function
    GridFunction<TDim, ControlPointType>& ControlPointGridFunction() {return *(this->pGetGridFunction(CONTROL_POINT));}

    /// Get the control point grid function
    const GridFunction<TDim, ControlPointType>& ControlPointGridFunction() const {return *(this->pGetGridFunction(CONTROL_POINT));}

    /// Get the control point grid function pointer
    typename GridFunction<TDim, ControlPointType>::Pointer pControlPointGridFunction() {return this->pGetGridFunction(CONTROL_POINT);}

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::ConstPointer pControlPointGridFunction() const {return this->pGetGridFunction(CONTROL_POINT);}

    /// Get the control point weights vector
    std::vector<double> GetControlWeights() const
    {
        typename ControlGrid<ControlPointType>::ConstPointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        std::vector<double> Weights(pControlPointGrid->size());
        for (std::size_t i = 0; i < pControlPointGrid->size(); ++i)
            Weights[i] = (*pControlPointGrid)[i].W();
        return Weights;
    }

    /// Apply the homogeneous transformation to the patch by applying the homogeneous transformation to the control points grid. For DISPLACEMENT, access the grid function for DISPLACEMENT directly and transform it.
    void ApplyTransformation(const TransformationType& trans)
    {
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        ControlGridUtility::ApplyTransformation(*pControlPointGrid, trans);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function
    template<typename TDataType>
    typename GridFunction<TDim, TDataType>::Pointer CreateGridFunction(typename ControlGrid<TDataType>::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename GridFunction<TDim, TDataType>::Pointer pNewGridFunc = GridFunction<TDim, TDataType>::Create(mFESpace, pControlGrid);
        mpGridFunctions.push_back(pNewGridFunc);
        return pNewGridFunc;
    }

    /// Create and add the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::Pointer CreateGridFunction(const TVariableType& rVariable,
            typename ControlGrid<typename TVariableType::Type>::Pointer pControlGrid)
    {
        return this->CreateGridFunction<typename TVariableType::Type>(pControlGrid);
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::Pointer pGetGridFunction(const TVariableType& rVariable)
    {
        typedef typename GridFunction<TDim, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (std::size_t i = 0; i < mpGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(mpGridFunctions[i]);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                    return pGridFunc;
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        std::stringstream ss;
        ss << "The grid function with control grid " << rVariable.Name() << " does not exist in the database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::ConstPointer pGetGridFunction(const TVariableType& rVariable) const
    {
        typedef typename GridFunction<TDim, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (std::size_t i = 0; i < mpGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(mpGridFunctions[i]);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                    return pGridFunc;
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        std::stringstream ss;
        ss << "The grid function with control grid " << rVariable.Name() << " does not exist in the database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    /// Filter out and get the underlying double grid functions
    DoubleGridFunctionContainerType DoubleGridFunctions() {return this->ExtractGridFunctions<DoubleGridFunctionContainerType, double>(mpGridFunctions);}
    DoubleGridFunctionContainerType DoubleGridFunctions() const {return this->ExtractGridFunctions<DoubleGridFunctionContainerType, double>(mpGridFunctions);}

    /// Filter out and get the underlying array_1d grid functions
    Array1DGridFunctionContainerType Array1DGridFunctions() {return this->ExtractGridFunctions<Array1DGridFunctionContainerType, array_1d<double, 3> >(mpGridFunctions);}
    Array1DGridFunctionContainerType Array1DGridFunctions() const {return this->ExtractGridFunctions<Array1DGridFunctionContainerType, array_1d<double, 3> >(mpGridFunctions);}

    /// Filter out and get the underlying Vector grid functions
    VectorGridFunctionContainerType VectorGridFunctions() {return this->ExtractGridFunctions<VectorGridFunctionContainerType, Vector>(mpGridFunctions);}
    VectorGridFunctionContainerType VectorGridFunctions() const {return this->ExtractGridFunctions<VectorGridFunctionContainerType, Vector>(mpGridFunctions);}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Validate the patch
    virtual bool Validate() const
    {
        if (Id() == 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The patch must have an Id", "")
        }

        if (pControlPointGridFunction() != NULL)
            if (pControlPointGridFunction()->pControlGrid()->Size() != this->TotalNumber())
                KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions().begin();
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
        for (int i = _LEFT_; i <= _BACK_; ++i)
        {
            BoundarySide side = static_cast<BoundarySide>(i);

            if (this->pNeighbor(side) != NULL)
            {
                if (this->Type() != this->pNeighbor(side)->Type())
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The patch type between this and the neighbor is incompatible", "")
                }
                else
                {
                    // find the side of the other neighbor
                    BoundarySide other_side = this->pNeighbor(side)->FindBoundarySide(this->shared_from_this());

                    if (other_side == _NUMBER_OF_BOUNDARY_SIDE)
                        KRATOS_THROW_ERROR(std::logic_error, "No neighbor of the neighbor is the same as this. Error setting the neighbor.", "")

                    bool check = CheckBoundaryCompatibility(*this, side, *(this->pNeighbor(side)), other_side);
                    if (!check)
                    {
                        KRATOS_WATCH(side)
                        KRATOS_WATCH(other_side)
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary between this and the neighbor is incompatible", "")
                        return false;
                    }
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

    /// Get/Set the neighbor
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

    /// Find the boundary side of the neighbor
    BoundarySide FindBoundarySide(typename Patch<TDim>::ConstPointer pPatch) const
    {
        for (int i = _LEFT_; i <= _BACK_; ++i)
        {
            BoundarySide side = static_cast<BoundarySide>(i);
            if (this->pNeighbor(side) == pPatch)
                return side;
        }
        return _NUMBER_OF_BOUNDARY_SIDE;
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
        rOStream << Type() << ", Id = " << Id() << ", Addr = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        if (pFESpace() != NULL)
            rOStream << *pFESpace() << std::endl;
        if (pControlPointGridFunction() != NULL)
            rOStream << *(pControlPointGridFunction()->pControlGrid());
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

    // container to contain all the grid functions
    std::vector<boost::any> mpGridFunctions; // using boost::any so store pointers to grid function

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

    /// Helper to extract the grid functions out from boost::any
    template<class TContainerType, typename TDataType>
    TContainerType ExtractGridFunctions(const std::vector<boost::any>& pGridFunctions) const
    {
        TContainerType GridFuncs;

        typedef typename GridFunction<TDim, TDataType>::Pointer GridFunctionPointerType;

        for (std::size_t i = 0; i < pGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(mpGridFunctions[i]);
                GridFuncs.push_back(pGridFunc);
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }

        return GridFuncs;
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
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
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
    MultiPatch() : mIsEnumerated(false) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPatch(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
        mIsEnumerated = false;
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

    /// Check the enumeration flag
    const bool& IsEnumerated() const {return mIsEnumerated;}

    /// Get the equation system size
    std::size_t EquationSystemSize() const {return mEquationSystemSize;}

    /// Locate the patch of the global equation id and the corresponding local id to determine the control value
    std::tuple<std::size_t, std::size_t> EquationIdLocation(const std::size_t& global_id) const
    {
        assert(IsEnumerated());

        std::map<std::size_t, std::size_t>::const_iterator it = mGlobalToPatch.find(global_id);
        if (it == mGlobalToPatch.end())
            KRATOS_THROW_ERROR(std::logic_error, "The global id does not exist in the global_to_patch map.", "")

        const std::size_t& patch_id = it->second;
        const std::size_t& local_id = pGetPatch(it->second)->pFESpace()->LocalId(global_id);

        return std::make_tuple(patch_id, local_id);
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

    /// Get the patch with specific Id
    typename PatchType::Pointer pGetPatch(const std::size_t& Id)
    {
        typename PatchContainerType::ptr_iterator it_patch = mpPatches.find(Id).base();
        assert(it_patch != mpPatches.ptr_end());
        return *it_patch;
    }

    /// Get the patch with specific Id
    typename PatchType::ConstPointer pGetPatch(const std::size_t& Id) const
    {
        typename PatchContainerType::ptr_const_iterator it_patch = mpPatches.find(Id).base();
        assert(it_patch != mpPatches.ptr_end());
        return typename PatchType::ConstPointer(*it_patch);
    }

    /// Access the underlying list of patches
    /// WARNING!!! be careful with this routine
    PatchContainerType& Patches() {return mpPatches;}

    /// Access the underlying list of patches
    const PatchContainerType& Patches() const {return mpPatches;}

    /// iterators
    typename PatchContainerType::iterator begin() {return mpPatches.begin();}
    typename PatchContainerType::const_iterator begin() const {return mpPatches.begin();}
    typename PatchContainerType::iterator end() {return mpPatches.end();}
    typename PatchContainerType::const_iterator end() const {return mpPatches.end();}

    /// Get the number of patches
    const std::size_t& size() const {return mpPatches.size();}

    /// Enumerate all the patches
    std::size_t Enumerate()
    {
        // reset global ids for each patch
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->ResetFunctionIndices();
        }

        // enumerate each patch
        mEquationSystemSize = 0;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            mEquationSystemSize = (*it)->pFESpace()->Enumerate(mEquationSystemSize);
            KRATOS_WATCH(mEquationSystemSize)

            // transfer the enumeration to neighbor boundary
            for (int i = _LEFT_; i <= _BACK_; ++i)
            {
                BoundarySide side = static_cast<BoundarySide>(i);

                if ((*it)->pNeighbor(side) != NULL)
                {
                    // find the side of the other neighbor
                    BoundarySide other_side = (*it)->pNeighbor(side)->FindBoundarySide(*it);

                    if (other_side == _NUMBER_OF_BOUNDARY_SIDE)
                        KRATOS_THROW_ERROR(std::logic_error, "No neighbor of the neighbor is the same as this. Error setting the neighbor.", "")

                    // check the boundary compatibility again
                    if (!(*it)->CheckBoundaryCompatibility(side, *((*it)->pNeighbor(side)), other_side))
                    {
                        KRATOS_WATCH(side)
                        KRATOS_WATCH(other_side)
                        KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with the neighbor is not satisfied", "")
                    }
                    else
                    {
                        std::vector<std::size_t> func_indices = (*it)->pFESpace()->ExtractBoundaryFunctionIndices(side);
                        (*it)->pNeighbor(side)->pFESpace()->AssignBoundaryFunctionIndices(other_side, func_indices);
                    }
                }
            }
        }

        // collect all the enumerated numbers and reassign with new to make it consecutive
        std::set<std::size_t> all_indices;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            all_indices.insert((*it)->pFESpace()->FunctionIndices().begin(), (*it)->pFESpace()->FunctionIndices().end());
        }

        std::map<std::size_t, std::size_t> new_indices;
        std::size_t cnt = 0;
        for (std::set<std::size_t>::iterator it = all_indices.begin(); it != all_indices.end(); ++it)
        {
            new_indices[*it] = cnt++;
        }

        // reassign the new indices to each patch
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->UpdateFunctionIndices(new_indices);
        }

        // rebuild the global to patch map
        mGlobalToPatch.clear();
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            const std::vector<std::size_t>& global_indices = (*it)->pFESpace()->FunctionIndices();
            for (std::size_t i = 0; i < global_indices.size(); ++i)
                mGlobalToPatch[global_indices[i]] = (*it)->Id();
        }

        // turn on the enumerated flag
        mIsEnumerated = true;

        return mEquationSystemSize;
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
    bool mIsEnumerated;
    std::size_t mEquationSystemSize;
    std::map<std::size_t, std::size_t> mGlobalToPatch; // this is to map each global id to a patch id

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


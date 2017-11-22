//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
// #include "custom_utilities/hierarchical_bsplines/hb_mesh.h"

namespace Kratos
{

/**
This class represents a single NURBS patch in parametric coordinates.
 */
class MultiPatchRefinementUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchRefinementUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiPatchRefinementUtility() {}

    /// Destructor
    virtual ~MultiPatchRefinementUtility() {}

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    template<int TDim>
    void InsertKnots(typename Patch<TDim>::Pointer& pPatch, const std::vector<std::vector<double> >& ins_knots)
    {
        std::set<std::size_t> refined_patches;
        InsertKnots<TDim>(pPatch, refined_patches, ins_knots);
    }

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    template<int TDim>
    void InsertKnots(typename Patch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::vector<double> >& ins_knots)
    {
        if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the NURBS patch")

        if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
        {
            // get the parent multipatch
            typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

            // create new patch with same Id
            typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

            // compute the transformation matrix
            Matrix T;
            std::vector<std::vector<double> > new_knots(TDim);
            std::vector<double> new_weights;
            std::vector<double> weights = pPatch->GetControlWeights();

            typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
            typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

            std::vector<std::size_t> new_size(TDim);

            this->ComputeNURBSKnotInsertionCoefficients<TDim>(T, new_knots, new_weights, pFESpace, ins_knots, weights);

            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                new_size[dim] = new_knots[dim].size() - pPatch->Order(dim) - 1;
                pNewFESpace->SetKnotVector(dim, new_knots[dim]);
                pNewFESpace->SetInfo(dim, new_size[dim], pPatch->Order(dim));
            }

            // KRATOS_WATCH(T)

            // set the new FESpace
            pNewPatch->SetFESpace(pNewFESpace);

            // transform and transfer the control points
            typename ControlGrid<ControlPoint<double> >::Pointer pNewControlPoints = typename ControlGrid<ControlPoint<double> >::Pointer (new RegularControlGrid<TDim, ControlPoint<double> >(new_size));
            ControlGridUtility::Transform<ControlPoint<double>, Matrix>(T, *(pPatch->pControlPointGridFunction()->pControlGrid()), *pNewControlPoints);
            pNewControlPoints->SetName(pPatch->pControlPointGridFunction()->pControlGrid()->Name());
            pNewPatch->CreateControlPointGridFunction(pNewControlPoints);

            // transfer the grid function
            for (typename Patch<TDim>::DoubleGridFunctionContainterType::const_iterator it = pPatch->DoubleGridFunctions().begin();
                    it != pPatch->DoubleGridFunctions().end(); ++it)
            {
                typename ControlGrid<double>::Pointer pNewDoubleControlGrid = typename ControlGrid<double>::Pointer (new RegularControlGrid<TDim, double>(new_size));
                ControlGridUtility::Transform<double, Matrix>(T, *((*it)->pControlGrid()), *pNewDoubleControlGrid);
                pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
                pNewPatch->CreateDoubleGridFunction(pNewDoubleControlGrid);
            }

            for (typename Patch<TDim>::Array1DGridFunctionContainerType::const_iterator it = pPatch->Array1DGridFunctions().begin();
                    it != pPatch->Array1DGridFunctions().end(); ++it)
            {
                typename ControlGrid<array_1d<double, 3> >::Pointer pNewArray1DControlGrid = typename ControlGrid<array_1d<double, 3> >::Pointer (new RegularControlGrid<TDim, array_1d<double, 3> >(new_size));
                ControlGridUtility::Transform<array_1d<double, 3>, Matrix>(T, *((*it)->pControlGrid()), *pNewArray1DControlGrid);
                pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
                pNewPatch->CreateArray1DGridFunction(pNewArray1DControlGrid);
            }

            for (typename Patch<TDim>::VectorGridFunctionContainerType::const_iterator it = pPatch->VectorGridFunctions().begin();
                    it != pPatch->VectorGridFunctions().end(); ++it)
            {
                typename ControlGrid<Vector>::Pointer pNewVectorControlGrid = typename ControlGrid<Vector>::Pointer (new RegularControlGrid<TDim, Vector>(new_size));
                ControlGridUtility::Transform<Vector, Matrix>(T, *((*it)->pControlGrid()), *pNewVectorControlGrid);
                pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
                pNewPatch->CreateVectorGridFunction(pNewVectorControlGrid);
            }

            // mark refined patch
            refined_patches.insert(pPatch->Id());

            // transfer the inserted knots to neighbors
            std::vector<std::vector<double> > neib_ins_knots(TDim);

            if (pPatch->pNeighbor(_LEFT_) != NULL)
            {
                if (pPatch->pNeighbor(_LEFT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_LEFT_, pPatch->pNeighbor(_LEFT_));
                    pPatch->pNeighbor(_LEFT_)->pSetNeighbor(_RIGHT_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        neib_ins_knots[2] = ins_knots[2];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_LEFT_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_LEFT_)->pFESpace()->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_RIGHT_) != NULL)
            {
                if (pPatch->pNeighbor(_RIGHT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_RIGHT_, pPatch->pNeighbor(_RIGHT_));
                    pPatch->pNeighbor(_RIGHT_)->pSetNeighbor(_LEFT_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        neib_ins_knots[2] = ins_knots[2];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_RIGHT_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_RIGHT_)->pFESpace()->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_TOP_) != NULL)
            {
                if (pPatch->pNeighbor(_TOP_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_TOP_, pPatch->pNeighbor(_TOP_));
                    pPatch->pNeighbor(_TOP_)->pSetNeighbor(_BOTTOM_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        neib_ins_knots[1] = ins_knots[1];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_TOP_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_TOP_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_BOTTOM_) != NULL)
            {
                if (pPatch->pNeighbor(_BOTTOM_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_BOTTOM_, pPatch->pNeighbor(_BOTTOM_));
                    pPatch->pNeighbor(_BOTTOM_)->pSetNeighbor(_TOP_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        neib_ins_knots[1] = ins_knots[1];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_BOTTOM_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_BOTTOM_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_FRONT_) != NULL)
            {
                if (pPatch->pNeighbor(_FRONT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_FRONT_, pPatch->pNeighbor(_FRONT_));
                    pPatch->pNeighbor(_FRONT_)->pSetNeighbor(_BACK_, pNewPatch);

                    neib_ins_knots[0] = ins_knots[0];
                    neib_ins_knots[2] = ins_knots[2];
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_FRONT_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_FRONT_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_BACK_) != NULL)
            {
                if (pPatch->pNeighbor(_BACK_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_BACK_, pPatch->pNeighbor(_BACK_));
                    pPatch->pNeighbor(_BACK_)->pSetNeighbor(_FRONT_, pNewPatch);

                    neib_ins_knots[0] = ins_knots[0];
                    neib_ins_knots[2] = ins_knots[2];
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_BACK_);
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                // else if (pPatch->pNeighbor(_BACK_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            // swap
            // KRATOS_WATCH(*pNewPatch)
            // std::cout << pPatch << " will be swapped with " << pNewPatch << std::endl;
            pPatch.swap(pNewPatch);
            // std::cout << "  Input patch now is " << pPatch << ", and the new patch becomes " << pNewPatch << std::endl;
            // KRATOS_WATCH(*pPatch)

            // replace the corresponding patch in multipatch
            pMultiPatch->Patches()(pPatch->Id()) = pPatch;

            // TODO re-generate the connection topology
            // 
        }
    }

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<int TDim>
    void DegreeElevate(typename Patch<TDim>::Pointer& pPatch, const std::vector<std::size_t>& order_increment)
    {
        std::set<std::size_t> refined_patches;
        DegreeElevate<TDim>(pPatch, refined_patches, order_increment);
    }

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<int TDim>
    void DegreeElevate(typename Patch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::size_t>& order_increment)
    {
        if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
        {
            // get the parent multipatch
            typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

            // create new patch with same Id
            typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

            // elevate the degree and initialize new patch
            typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
            typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

            std::vector<std::vector<double> > new_knots(TDim);

            std::vector<std::size_t> new_size(TDim);
            for (std::size_t i = 0; i < TDim; ++i)
                new_size[i] = pFESpace->Number(i);

            typename RegularControlGrid<TDim, ControlPoint<double> >::Pointer pControlPoints
                = boost::dynamic_pointer_cast<RegularControlGrid<TDim, ControlPoint<double> > >(pPatch->pControlPointGridFunction()->pControlGrid());

            typename RegularControlGrid<TDim, ControlPoint<double> >::Pointer pNewControlPoints
                = typename RegularControlGrid<TDim, ControlPoint<double> >::Pointer(new RegularControlGrid<TDim, ControlPoint<double> >(new_size)); // note here that the size is just temporary, it will be raised later on.

            this->ComputeBsplinesDegreeElevation<TDim>(*pControlPoints, *pFESpace, order_increment, *pNewControlPoints, new_knots);

            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                new_size[dim] = new_knots[dim].size() - pFESpace->Order(dim) - order_increment[dim] - 1;
                pNewFESpace->SetKnotVector(dim, new_knots[dim]);
                pNewFESpace->SetInfo(dim, new_size[dim], pFESpace->Order(dim) + order_increment[dim]);
            }

            pNewControlPoints->SetName(pPatch->pControlPointGridFunction()->pControlGrid()->Name());
            pNewPatch->SetFESpace(pNewFESpace);
            pNewPatch->CreateControlPointGridFunction(pNewControlPoints);

            // mark refined patch
            refined_patches.insert(pPatch->Id());

            // transfer the order increment to neighbors
            std::vector<std::size_t> neib_order_increment(TDim);

            if (pPatch->pNeighbor(_LEFT_) != NULL)
            {
                if (pPatch->pNeighbor(_LEFT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_LEFT_, pPatch->pNeighbor(_LEFT_));
                    pPatch->pNeighbor(_LEFT_)->pSetNeighbor(_RIGHT_, pNewPatch);

                    if(TDim == 1)
                    {
                        neib_order_increment[0] = order_increment[0];
                    }
                    else if(TDim == 2)
                    {
                        neib_order_increment[1] = order_increment[1];
                    }
                    else if (TDim == 3)
                    {
                        neib_order_increment[1] = order_increment[1];
                        neib_order_increment[2] = order_increment[2];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_LEFT_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_LEFT_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_RIGHT_) != NULL)
            {
                if (pPatch->pNeighbor(_RIGHT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_RIGHT_, pPatch->pNeighbor(_RIGHT_));
                    pPatch->pNeighbor(_RIGHT_)->pSetNeighbor(_LEFT_, pNewPatch);

                    if(TDim == 1)
                    {
                        neib_order_increment[0] = order_increment[0];
                    }
                    else if (TDim == 2)
                    {
                        neib_order_increment[1] = order_increment[1];
                    }
                    else if (TDim == 3)
                    {
                        neib_order_increment[1] = order_increment[1];
                        neib_order_increment[2] = order_increment[2];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_RIGHT_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_RIGHT_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_TOP_) != NULL)
            {
                if (pPatch->pNeighbor(_TOP_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_TOP_, pPatch->pNeighbor(_TOP_));
                    pPatch->pNeighbor(_TOP_)->pSetNeighbor(_BOTTOM_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_order_increment[0] = order_increment[0];
                    }
                    else if (TDim == 3)
                    {
                        neib_order_increment[0] = order_increment[0];
                        neib_order_increment[1] = order_increment[1];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_TOP_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_TOP_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_BOTTOM_) != NULL)
            {
                if (pPatch->pNeighbor(_BOTTOM_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_BOTTOM_, pPatch->pNeighbor(_BOTTOM_));
                    pPatch->pNeighbor(_BOTTOM_)->pSetNeighbor(_TOP_, pNewPatch);

                    if (TDim == 2)
                    {
                        neib_order_increment[0] = order_increment[0];
                    }
                    else if (TDim == 3)
                    {
                        neib_order_increment[0] = order_increment[0];
                        neib_order_increment[1] = order_increment[1];
                    }
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_BOTTOM_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_BOTTOM_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_FRONT_) != NULL)
            {
                if (pPatch->pNeighbor(_FRONT_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_FRONT_, pPatch->pNeighbor(_FRONT_));
                    pPatch->pNeighbor(_FRONT_)->pSetNeighbor(_BACK_, pNewPatch);

                    neib_order_increment[0] = order_increment[0];
                    neib_order_increment[2] = order_increment[2];
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_FRONT_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_FRONT_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            if (pPatch->pNeighbor(_BACK_) != NULL)
            {
                if (pPatch->pNeighbor(_BACK_)->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
                {
                    pNewPatch->pSetNeighbor(_BACK_, pPatch->pNeighbor(_BACK_));
                    pPatch->pNeighbor(_BACK_)->pSetNeighbor(_FRONT_, pNewPatch);

                    neib_order_increment[0] = order_increment[0];
                    neib_order_increment[2] = order_increment[2];
                    typename Patch<TDim>::Pointer pNeighbor = pPatch->pNeighbor(_BACK_);
                    DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);
                }
                // else if (pPatch->pNeighbor(_BACK_)->Type() == HBMesh<TDim>::StaticType())
                // {
                //     //TODO
                //     KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                // }
            }

            // swap
            // KRATOS_WATCH(*pNewPatch)
            // std::cout << pPatch << " will be swapped with " << pNewPatch << std::endl;
            pPatch.swap(pNewPatch);
            // std::cout << "  Input patch now is " << pPatch << ", and the new patch becomes " << pNewPatch << std::endl;
            // KRATOS_WATCH(*pPatch)

            // replace the corresponding patch in multipatch
            pMultiPatch->Patches()(pPatch->Id()) = pPatch;
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchRefinementUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    /// Compute the transformation matrix for knot insertion (NURBS version)
    template<int TDim>
    void ComputeNURBSKnotInsertionCoefficients(
        Matrix& T,
        std::vector<std::vector<double> >& new_knots,
        std::vector<double>& new_weights,
        typename BSplinesFESpace<TDim>::Pointer& pFESpace,
        const std::vector<std::vector<double> >& ins_knots,
        const std::vector<double>& weights) const
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not implemented for dimension " << TDim;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    template<int TDim>
    void ComputeBsplinesDegreeElevation(
        const RegularControlGrid<TDim, ControlPoint<double> >& ControlPoints,
        const BSplinesFESpace<TDim>& rFESpace,
        const std::vector<std::size_t>& order_increment,
        RegularControlGrid<TDim, ControlPoint<double> >& NewControlPoints,
        std::vector<std::vector<double> >& new_knots) const
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not implemented for dimension " << TDim;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

};

template<>
void MultiPatchRefinementUtility::ComputeNURBSKnotInsertionCoefficients<1>(
    Matrix& T,
    std::vector<std::vector<double> >& new_knots,
    std::vector<double>& new_weights,
    typename BSplinesFESpace<1>::Pointer& pFESpace,
    const std::vector<std::vector<double> >& ins_knots,
    const std::vector<double>& weights) const
{
    BSplineUtils::ComputeNURBSKnotInsertionCoefficients1D(T,
            new_knots[0],
            new_weights,
            pFESpace->Order(0),
            pFESpace->KnotVector(0),
            ins_knots[0],
            weights);
}

template<>
void MultiPatchRefinementUtility::ComputeNURBSKnotInsertionCoefficients<2>(
    Matrix& T,
    std::vector<std::vector<double> >& new_knots,
    std::vector<double>& new_weights,
    typename BSplinesFESpace<2>::Pointer& pFESpace,
    const std::vector<std::vector<double> >& ins_knots,
    const std::vector<double>& weights) const
{
    BSplineUtils::ComputeNURBSKnotInsertionCoefficients2D(T,
            new_knots[0], new_knots[1],
            new_weights,
            pFESpace->Order(0), pFESpace->Order(1),
            pFESpace->KnotVector(0), pFESpace->KnotVector(1),
            ins_knots[0], ins_knots[1],
            weights);
}

template<>
void MultiPatchRefinementUtility::ComputeNURBSKnotInsertionCoefficients<3>(
    Matrix& T,
    std::vector<std::vector<double> >& new_knots,
    std::vector<double>& new_weights,
    typename BSplinesFESpace<3>::Pointer& pFESpace,
    const std::vector<std::vector<double> >& ins_knots,
    const std::vector<double>& weights) const
{
    BSplineUtils::ComputeNURBSKnotInsertionCoefficients3D(T,
            new_knots[0], new_knots[1], new_knots[2],
            new_weights,
            pFESpace->Order(0), pFESpace->Order(1), pFESpace->Order(2),
            pFESpace->KnotVector(0), pFESpace->KnotVector(1), pFESpace->KnotVector(2),
            ins_knots[0], ins_knots[1], ins_knots[2],
            weights);
}

template<>
void MultiPatchRefinementUtility::ComputeBsplinesDegreeElevation<1>(
    const RegularControlGrid<1, ControlPoint<double> >& ControlPoints,
    const BSplinesFESpace<1>& rFESpace,
    const std::vector<std::size_t>& order_increment,
    RegularControlGrid<1, ControlPoint<double> >& NewControlPoints,
    std::vector<std::vector<double> >& new_knots) const
{
    ControlPoint<double> null_control_point(0.0);

    BSplineUtils::ComputeBsplinesDegreeElevation1D(rFESpace.Order(0),
            ControlPoints,
            rFESpace.KnotVector(0),
            order_increment[0],
            NewControlPoints,
            new_knots[0],
            null_control_point);
}

template<>
void MultiPatchRefinementUtility::ComputeBsplinesDegreeElevation<2>(
    const RegularControlGrid<2, ControlPoint<double> >& ControlPoints,
    const BSplinesFESpace<2>& rFESpace,
    const std::vector<std::size_t>& order_increment,
    RegularControlGrid<2, ControlPoint<double> >& NewControlPoints,
    std::vector<std::vector<double> >& new_knots) const
{
    ControlPoint<double> null_control_point(0.0);

    BSplineUtils::ComputeBsplinesDegreeElevation2D(rFESpace.Order(0), rFESpace.Order(1),
            ControlPoints,
            rFESpace.KnotVector(0), rFESpace.KnotVector(1),
            order_increment[0], order_increment[1],
            NewControlPoints,
            new_knots[0], new_knots[1],
            null_control_point);
}

template<>
void MultiPatchRefinementUtility::ComputeBsplinesDegreeElevation<3>(
    const RegularControlGrid<3, ControlPoint<double> >& ControlPoints,
    const BSplinesFESpace<3>& rFESpace,
    const std::vector<std::size_t>& order_increment,
    RegularControlGrid<3, ControlPoint<double> >& NewControlPoints,
    std::vector<std::vector<double> >& new_knots) const
{
    ControlPoint<double> null_control_point(0.0);

    BSplineUtils::ComputeBsplinesDegreeElevation3D(rFESpace.Order(0), rFESpace.Order(1), rFESpace.Order(2),
            ControlPoints,
            rFESpace.KnotVector(0), rFESpace.KnotVector(1), rFESpace.KnotVector(2),
            order_increment[0], order_increment[1], order_increment[2],
            NewControlPoints,
            new_knots[0], new_knots[1], new_knots[2],
            null_control_point);
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchRefinementUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED defined


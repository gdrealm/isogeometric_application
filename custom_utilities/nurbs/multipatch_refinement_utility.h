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
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/grid_function.h"
#include "custom_utilities/nurbs/grid_function_utility.h"
#include "custom_utilities/nurbs/patch.h"
#include "custom_utilities/nurbs/nurbs_patch.h"
#include "custom_utilities/hierarchical_nurbs/hn_mesh.h"

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
    template<std::size_t TDim>
    static void InsertKnots(typename NURBSPatch<TDim>::Pointer& pPatch, const std::vector<std::vector<double> >& ins_knots)
    {
        std::set<std::size_t> refined_patches;
        InsertKnots<TDim>(pPatch, refined_patches, ins_knots);
    }

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    template<std::size_t TDim>
    static void InsertKnots(typename NURBSPatch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::vector<double> >& ins_knots)
    {
        if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
        {
            // get the parent multipatch
            typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

            // create new patch with same Id
            typename NURBSPatch<TDim>::Pointer pNewPatch = typename NURBSPatch<TDim>::Pointer(new NURBSPatch<TDim>(pPatch->Id()));

            // compute the transformation matrix
            Matrix T;
            std::vector<std::vector<double> > new_knots(TDim);
            std::vector<double> new_weights;
            std::vector<double> weights = pPatch->GetControlWeights();

            std::vector<std::size_t> new_size(TDim);

            if (TDim == 1)
            {
                BSplineUtils::ComputeNURBSKnotInsertionCoefficients1D(T, new_knots[0], new_weights,
                        pPatch->Order(0),
                        pPatch->KnotVector(0), ins_knots[0], weights);
            }
            else if (TDim == 2)
            {
                BSplineUtils::ComputeNURBSKnotInsertionCoefficients2D(T, new_knots[0], new_knots[1], new_weights,
                        pPatch->Order(0), pPatch->Order(1),
                        pPatch->KnotVector(0), pPatch->KnotVector(1), ins_knots[0], ins_knots[1], weights);
            }
            else if (TDim == 3)
            {
                BSplineUtils::ComputeNURBSKnotInsertionCoefficients3D(T, new_knots[0], new_knots[1], new_knots[2], new_weights,
                        pPatch->Order(0), pPatch->Order(1), pPatch->Order(2),
                        pPatch->KnotVector(0), pPatch->KnotVector(1), pPatch->KnotVector(2), ins_knots[0], ins_knots[1], ins_knots[2], weights);
            }

            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                new_size[dim] = new_knots[dim].size() - pPatch->Order(dim) - 1;
                pNewPatch->SetKnotVector(dim, new_knots[dim]);
                pNewPatch->SetInfo(dim, new_size[dim], pPatch->Order(dim));
            }

            // KRATOS_WATCH(T)

            // transform and transfer the control points
            typename GridFunction<TDim, ControlPoint<double> >::Pointer pNewControlPoints
                = GridFunctionUtility::Transform<TDim, ControlPoint<double>, Matrix>(T, new_size, pPatch->ControlPointGrid());
            pNewControlPoints->SetName(pPatch->ControlPointGrid()->Name());
            pNewPatch->SetControlPointGrid(pNewControlPoints);

            // transfer the grid function
            for (typename Patch<TDim>::DoubleGridFunctionContainterType::const_iterator it = pPatch->DoubleGridFunctions().begin();
                    it != pPatch->DoubleGridFunctions().end(); ++it)
            {
                typename GridFunction<TDim, double>::Pointer pNewDoubleGridFunction
                    = GridFunctionUtility::Transform<TDim, double, Matrix>(T, new_size, *it);
                pNewDoubleGridFunction->SetName((*it)->Name());
                pNewPatch->AddDoubleGridFunction(pNewDoubleGridFunction);
            }

            for (typename Patch<TDim>::Array1DGridFunctionContainerType::const_iterator it = pPatch->Array1DGridFunctions().begin();
                    it != pPatch->Array1DGridFunctions().end(); ++it)
            {
                typename GridFunction<TDim, array_1d<double, 3> >::Pointer pNewArray1DGridFunction
                    = GridFunctionUtility::Transform<TDim, array_1d<double, 3>, Matrix>(T, new_size, *it);
                pNewArray1DGridFunction->SetName((*it)->Name());
                pNewPatch->AddArray1DGridFunction(pNewArray1DGridFunction);
            }

            for (typename Patch<TDim>::VectorGridFunctionContainerType::const_iterator it = pPatch->VectorGridFunctions().begin();
                    it != pPatch->VectorGridFunctions().end(); ++it)
            {
                typename GridFunction<TDim, Vector>::Pointer pNewVectorGridFunction
                    = GridFunctionUtility::Transform<TDim, Vector, Matrix>(T, new_size, *it);
                pNewVectorGridFunction->SetName((*it)->Name());
                pNewPatch->AddVectorGridFunction(pNewVectorGridFunction);
            }

            // mark refined patch
            refined_patches.insert(pPatch->Id());

            // transfer the inserted knots to neighbors
            std::vector<std::vector<double> > neib_ins_knots(TDim);

            if (pPatch->pLeft() != NULL)
            {
                if (pPatch->pLeft()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetLeft(pPatch->pLeft());
                    pPatch->pLeft()->pSetRight(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pLeft());
                    if (TDim == 2)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        neib_ins_knots[2] = ins_knots[2];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                }
                else if (pPatch->pLeft()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
            }

            if (pPatch->pRight() != NULL)
            {
                if (pPatch->pRight()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetRight(pPatch->pRight());
                    pPatch->pRight()->pSetLeft(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pRight());

                    if (TDim == 2)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[1] = ins_knots[1];
                        neib_ins_knots[2] = ins_knots[2];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                }
                else if (pPatch->pRight()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
            }

            if (pPatch->pTop() != NULL)
            {
                if (pPatch->pTop()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetTop(pPatch->pTop());
                    pPatch->pTop()->pSetBottom(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pTop());
                    if (TDim == 2)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        neib_ins_knots[1] = ins_knots[1];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                }
                else if (pPatch->pTop()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
            }

            if (pPatch->pBottom() != NULL)
            {
                if (pPatch->pBottom()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetBottom(pPatch->pBottom());
                    pPatch->pBottom()->pSetTop(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pBottom());
                    if (TDim == 2)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                    else if (TDim == 3)
                    {
                        neib_ins_knots[0] = ins_knots[0];
                        neib_ins_knots[1] = ins_knots[1];
                        InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                    }
                }
                else if (pPatch->pBottom()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
            }

            if (pPatch->pFront() != NULL)
            {
                if (pPatch->pFront()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetFront(pPatch->pFront());
                    pPatch->pFront()->pSetBack(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pFront());
                    neib_ins_knots[0] = ins_knots[0];
                    neib_ins_knots[2] = ins_knots[2];
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                else if (pPatch->pFront()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
            }

            if (pPatch->pBack() != NULL)
            {
                if (pPatch->pBack()->Type() == NURBSPatch<TDim>::StaticType())
                {
                    pNewPatch->pSetBack(pPatch->pBack());
                    pPatch->pBack()->pSetFront(pNewPatch);

                    typename NURBSPatch<TDim>::Pointer pNeighbor = boost::dynamic_pointer_cast<NURBSPatch<TDim> >(pPatch->pBack());
                    neib_ins_knots[0] = ins_knots[0];
                    neib_ins_knots[2] = ins_knots[2];
                    InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);
                }
                else if (pPatch->pBack()->Type() == HnMesh<TDim>::StaticType())
                {
                    //TODO
                    KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
                }
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

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<std::size_t TDim>
    static void DegreeElevate(typename NURBSPatch<TDim>::Pointer& pPatch, const std::vector<std::size_t>& order_increment)
    {
        std::set<std::size_t> refined_patches;
        DegreeElevate<TDim>(pPatch, refined_patches, order_increment);
    }

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<std::size_t TDim>
    static void DegreeElevate(typename NURBSPatch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::size_t>& order_increment)
    {
        if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
        {
            // get the parent multipatch
            typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

            // create new patch with same Id
            typename NURBSPatch<TDim>::Pointer pNewPatch = typename NURBSPatch<TDim>::Pointer(new NURBSPatch<TDim>(pPatch->Id()));

            std::vector<std::vector<double> > new_knots(TDim);
            std::vector<double> new_weights;
            std::vector<double> weights = pPatch->GetControlWeights();

            std::vector<std::size_t> new_size(TDim);

            typename GridFunction<TDim, ControlPoint<double> >::Pointer pNewControlPoints
                = typename GridFunction<TDim, ControlPoint<double> >::Pointer(new GridFunction<TDim, ControlPoint<double> >(1));

            ControlPoint<double> null_control_point(0.0);

            if (TDim == 1)
            {
                KRATOS_WATCH(order_increment[0])
                BSplineUtils::ComputeBsplinesDegreeElevation1D(pPatch->Order(0),
                        *(pPatch->ControlPointGrid()),
                        pPatch->KnotVector(0),
                        order_increment[0],
                        *pNewControlPoints,
                        new_knots[0],
                        null_control_point);
                KRATOS_WATCH(__LINE__)
            }

            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                new_size[dim] = new_knots[dim].size() - pPatch->Order(dim) - 1;
                pNewPatch->SetKnotVector(dim, new_knots[dim]);
                pNewPatch->SetInfo(dim, new_size[dim], pPatch->Order(dim));
            }

            pNewControlPoints->SetName(pPatch->ControlPointGrid()->Name());
            pNewPatch->SetControlPointGrid(pNewControlPoints);

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

};

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


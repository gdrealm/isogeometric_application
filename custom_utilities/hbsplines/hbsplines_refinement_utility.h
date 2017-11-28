//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED

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
#include "custom_utilities/hbsplines/hbsplines_fespace.h"

namespace Kratos
{

/**
Class accounts for linear independent refinement of a single hierarchical B-Splines patch
 */
class HBSplinesRefinementUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesRefinementUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    HBSplinesRefinementUtility() {}

    /// Destructor
    virtual ~HBSplinesRefinementUtility() {}

    /// Refine a single B-Splines basis function
    template<int TDim>
    static void Refine(typename Patch<TDim>::Pointer pPatch, const std::size_t& Id, const int& EchoLevel)
    {
        if (pPatch->pFESpace()->Type() != HBSplinesFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

        // Type definitions
        typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;
        typedef typename HBSplinesFESpace<TDim>::bf_container_t bf_container_t;
        typedef typename HBSplinesFESpace<TDim>::CellType CellType;
        typedef typename HBSplinesFESpace<TDim>::cell_t cell_t;
        typedef typename HBSplinesFESpace<TDim>::cell_container_t cell_container_t;
        typedef typename Patch<TDim>::ControlPointType ControlPointType;

        // TODO get the list of variables in the patch so that it can be transferred during refinement


        // extract the hierarchical B-Splines space
        typename HBSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pPatch->pFESpace());

        #ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
        #endif

        bf_t p_bf;
        bool found = false;
        for(typename bf_container_t::iterator it = pFESpace->bf_begin(); it != pFESpace->bf_end(); ++it)
        {
            if((*it)->Id() == Id)
            {
                p_bf = *it;
                found = true;
            }
        }

        if(!found) return;

        // does not refine if maximum level is reached
        if(p_bf->Level() == pFESpace->MaxLevel())
        {
            std::cout << "Maximum level is reached, basis function " << p_bf->Id() << " is skipped" << std::endl;
            return;
        }

        /* create a list of basis function in the next level representing this basis function */
        // create a list of new knots
        double tol = 1.0e-10; // TODO what is this? can we parameterize?
        std::vector<std::vector<knot_t> > pnew_local_knots(TDim);
        std::vector<std::vector<double> > ins_knots(TDim);
        for(unsigned int dim = 0; dim < TDim; ++dim)
        {
            const std::vector<knot_t>& pLocalKnots = p_bf->LocalKnots(dim);

            for(std::vector<knot_t>::const_iterator it = pLocalKnots.begin(); it != pLocalKnots.end(); ++it)
            {
                pnew_local_knots[dim].push_back(*it);

                std::vector<knot_t>::const_iterator it2 = it + 1;
                if(it2 != pLocalKnots.end())
                {
                    if(fabs((*it2)->Value() - (*it)->Value()) > tol)
                    {
                        // now we just add the middle one, but in the general we can add arbitrary values
                        // TODO: find the way to generalize this or parameterize this
                        double ins_knot = 0.5 * ((*it)->Value() + (*it2)->Value());
                        knot_t p_new_knot;
                        p_new_knot = pFESpace->KnotVector(dim).pCreateUniqueKnot(ins_knot, tol);
                        pnew_local_knots[dim].push_back(p_new_knot);
                        ins_knots[dim].push_back(p_new_knot->Value());
                    }
                }
            }
        }

        /* compute the refinement coefficients */
        Vector RefinedCoeffs;
        std::vector<std::vector<double> > local_knots(TDim);
        for(std::size_t dim = 0; dim < TDim; ++dim)
            p_bf->GetLocalKnots(dim, local_knots[dim]);

        std::vector<std::vector<double> > new_knots(TDim);
        if (TDim == 2)
        {
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients2DLocal(RefinedCoeffs,
                new_knots[0], new_knots[1],
                pFESpace->Order(0), pFESpace->Order(1),
                local_knots[0], local_knots[1],
                ins_knots[0], ins_knots[1]);
        }
        else if (TDim == 3)
        {
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients3DLocal(RefinedCoeffs,
                new_knots[0], new_knots[1], new_knots[2],
                pFESpace->Order(0), pFESpace->Order(1), pFESpace->Order(2),
                local_knots[0], local_knots[1], local_knots[2],
                ins_knots[0], ins_knots[1], ins_knots[2]);
        }

        #ifdef ENABLE_PROFILING
        double time_1 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* create new basis functions */
        unsigned int next_level = p_bf->Level() + 1;
        if(next_level > pFESpace->LastLevel()) pFESpace->SetLastLevel(next_level);
        double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it.
        std::size_t last_id = pFESpace->LastId();
        typename cell_container_t::Pointer pnew_cells;

        if (TDim == 2)
        {
            pnew_cells = typename cell_container_t::Pointer(new CellManager2D<CellType>());

            std::size_t number_1 = pnew_local_knots[0].size() - pFESpace->Order(0) - 1;
            std::size_t number_2 = pnew_local_knots[1].size() - pFESpace->Order(1) - 1;
            for(std::size_t j = 0; j < number_2; ++j)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots2;
                for(std::size_t k = 0; k < pFESpace->Order(1) + 2; ++k)
                    pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

                for(std::size_t i = 0; i < number_1; ++i)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots1;
                    for(std::size_t k = 0; k < pFESpace->Order(0) + 2; ++k)
                        pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                    // create the basis function object
                    std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2};
                    bf_t pnew_bf = pFESpace->CreateBf(last_id+1, next_level, pLocalKnots);
                    if (pnew_bf->Id() == last_id+1) ++last_id;

                    std::size_t i_func = j * number_1 + i;

                    // transfer the control point information
                    ControlPointType oldC = p_bf->GetValue(CONTROL_POINT);
                    ControlPointType newC = pnew_bf->GetValue(CONTROL_POINT);
                    newC += RefinedCoeffs[i_func] * oldC;
                    pnew_bf->SetValue(CONTROL_POINT, newC);

                    // TODO transfer other control values from p_bf to pnew_bf

                    // create the cells for the basis function
                    for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
                    {
                        knot_t pLeft = pnew_local_knots[0][i + i1];
                        knot_t pRight = pnew_local_knots[0][i + i1 + 1];
                        for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                        {
                            knot_t pDown = pnew_local_knots[1][j + j1];
                            knot_t pUp = pnew_local_knots[1][j + j1 + 1];

                            // check if the cell domain area is nonzero
                            double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                            if(fabs(area) > area_tol)
                            {
                                std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                                cell_t pnew_cell = pFESpace->pCellManager()->CreateCell(pKnots);
                                pnew_cell->SetLevel(next_level);
                                pnew_bf->AddCell(pnew_cell);
                                pnew_cell->AddBf(pnew_bf);
                                pnew_cells->insert(pnew_cell);
                            }
                        }
                    }
                }
            }
        }
        else if (TDim == 3)
        {
            pnew_cells = typename cell_container_t::Pointer(new CellManager3D<CellType>());

            std::size_t number_1 = pnew_local_knots[0].size() - pFESpace->Order(0) - 1;
            std::size_t number_2 = pnew_local_knots[1].size() - pFESpace->Order(1) - 1;
            std::size_t number_3 = pnew_local_knots[2].size() - pFESpace->Order(2) - 1;
            for(std::size_t l = 0; l < number_3; ++l)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots3;
                for(std::size_t k = 0; k < pFESpace->Order(2) + 2; ++k)
                    pLocalKnots3.push_back(pnew_local_knots[2][l + k]);

                for(std::size_t j = 0; j < number_2; ++j)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots2;
                    for(std::size_t k = 0; k < pFESpace->Order(1) + 2; ++k)
                        pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

                    for(std::size_t i = 0; i < number_1; ++i)
                    {
                        // create and fill the local knot vector
                        std::vector<knot_t> pLocalKnots1;
                        for(std::size_t k = 0; k < pFESpace->Order(0) + 2; ++k)
                            pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                        // create the basis function object
                        std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2, pLocalKnots3};
                        bf_t pnew_bf = pFESpace->CreateBf(last_id+1, next_level, pLocalKnots);
                        if (pnew_bf->Id() == last_id+1) ++last_id;

                        // update the coordinates
                        std::size_t i_func = (l * number_2 + j) * number_1 + i;

                        // transfer the control point information
                        ControlPointType oldC = p_bf->GetValue(CONTROL_POINT);
                        ControlPointType newC = pnew_bf->GetValue(CONTROL_POINT);
                        newC += RefinedCoeffs[i_func] * oldC;
                        pnew_bf->SetValue(CONTROL_POINT, newC);

                        // TODO transfer other control values from p_bf to pnew_bf

                        // create the cells for the basis function
                        for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
                        {
                            knot_t pLeft = pnew_local_knots[0][i + i1];
                            knot_t pRight = pnew_local_knots[0][i + i1 + 1];

                            for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                            {
                                knot_t pDown = pnew_local_knots[1][j + j1];
                                knot_t pUp = pnew_local_knots[1][j + j1 + 1];

                                for(std::size_t l1 = 0; l1 < pFESpace->Order(2) + 1; ++l1)
                                {
                                    knot_t pBelow = pnew_local_knots[2][l + l1];
                                    knot_t pAbove = pnew_local_knots[2][l + l1 + 1];

                                    // check if the cell domain area is nonzero
                                    double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value()) * (pAbove->Value() - pBelow->Value());
                                    if(fabs(area) > area_tol)
                                    {
                                        std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp, pBelow, pAbove};
                                        cell_t pnew_cell = pFESpace->pCellManager()->CreateCell(pKnots);
                                        pnew_cell->SetLevel(next_level);
                                        pnew_bf->AddCell(pnew_cell);
                                        pnew_cell->AddBf(pnew_bf);
                                        pnew_cells->insert(pnew_cell);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        #ifdef ENABLE_PROFILING
        double time_2 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* remove the cells of the old basis function (remove only the cell in the current level) */
        typename cell_container_t::Pointer pcells_to_remove;
        if (TDim == 2)
            pcells_to_remove = typename cell_container_t::Pointer(new CellManager2D<CellType>());
        else if (TDim == 3)
            pcells_to_remove = typename cell_container_t::Pointer(new CellManager3D<CellType>());

        // firstly we check if the cell c of the current bf in the current level cover any sub-cells. Then the sub-cell includes all bfs of the cell c.
        for(typename HBSplinesFESpace<TDim>::BasisFunctionType::cell_iterator it_cell = p_bf->cell_begin(); it_cell != p_bf->cell_end(); ++it_cell)
        {
            if((*it_cell)->Level() == p_bf->Level())
            {
                for(typename cell_container_t::iterator it_subcell = pnew_cells->begin(); it_subcell != pnew_cells->end(); ++it_subcell)
                {
                    if((*it_subcell)->template IsCovered<TDim>(*it_cell))
                    {
                        for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                        {
                            (*it_subcell)->AddBf(*it_bf);
                            (*it_bf)->AddCell(*it_subcell);
                        }
                    }
                }

                // mark to remove the old cell
                pcells_to_remove->insert(*it_cell);
            }
        }

        // secondly, it happens that new cell c cover several existing cells. In this case cell c must be removed, and its bfs will be transferred to sub-cells.
        for(typename cell_container_t::iterator it_cell = pnew_cells->begin(); it_cell != pnew_cells->end(); ++it_cell)
        {
            std::vector<cell_t> p_cells = pFESpace->pCellManager()->GetCells(*it_cell);
            if(p_cells.size() > 0)
            {
                if((EchoLevel & ECHO_REFIMENT) == ECHO_REFIMENT)
                {
                    std::cout << "cell " << (*it_cell)->Id() << " is detected to contain some smaller cells:";
                    for(std::size_t i = 0; i < p_cells.size(); ++i)
                        std::cout << " " << p_cells[i]->Id();
                    std::cout << std::endl;
                }

                pcells_to_remove->insert(*it_cell);
                for(std::size_t i = 0; i < p_cells.size(); ++i)
                {
                    for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                    {
                        p_cells[i]->AddBf(*it_bf);
                        (*it_bf)->AddCell(p_cells[i]);
                    }
                }
            }
        }

        /* remove the cells */
        for(typename cell_container_t::iterator it_cell = pcells_to_remove->begin(); it_cell != pcells_to_remove->end(); ++it_cell)
        {
            pFESpace->pCellManager()->erase(*it_cell);
            for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
                (*it_bf)->RemoveCell(*it_cell);
        }

        /* remove the basis function from all the cells */
        for(typename cell_container_t::iterator it_cell = pFESpace->pCellManager()->begin(); it_cell != pFESpace->pCellManager()->end(); ++it_cell)
            (*it_cell)->RemoveBf(p_bf);

        /* remove the old basis function */
        pFESpace->RemoveBf(p_bf);

        #ifdef ENABLE_PROFILING
        double time_3 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* Debug part: to be removed when the code is stable
        for(typename cell_container_t::iterator it_cell = pFESpace->pCellManager()->begin(); it_cell != pFESpace->pCellManager()->end(); ++it_cell)
        {
            std::cout << "cell " << (*it_cell)->Id() << " supports:";
            for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                std::cout << " " << (*it_bf)->Id();
            std::cout << std::endl;
        }

        for(typename bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            std::cout << "bf " << (*it_bf)->Id() << " contains cell:";
            for(DeprecatedHBBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                std::cout << " " << (*it_cell)->Id();
            std::cout << std::endl;
        }
        */

        pFESpace->RecordRefinementHistory(p_bf->Id());
        if((EchoLevel & ECHO_REFIMENT) == ECHO_REFIMENT)
        {
            #ifdef ENABLE_PROFILING
            std::cout << "Refine bf " << p_bf->Id() << " completed" << std::endl;
            std::cout << " Time to compute the refinement coefficients: " << time_1 << " s" << std::endl;
            std::cout << " Time to create new cells and new bfs: " << time_2 << " s" << std::endl;
            std::cout << " Time to clean up: " << time_3 << " s" << std::endl;
            #else
            std::cout << "Refine bf " << p_bf->Id() << " completed" << std::endl;
            #endif
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HBSplinesRefinementUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesRefinementUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED defined


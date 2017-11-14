/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/nurbs/patch.h"
#include "custom_utilities/nurbs/nurbs_patch.h"
#include "custom_utilities/nurbs/nurbs_patch_library.h"
#include "custom_utilities/nurbs/nurbs_patch_utility.h"
#include "custom_utilities/nurbs/grid_function_utility.h"
#include "custom_utilities/nurbs/multipatch_refinement_utility.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_utilities/hierarchical_nurbs/hn_mesh.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

NURBSPatch<1>::Pointer NURBSPatchLibrary_CreateLinearPatch(NURBSPatchLibrary& rDummy, const std::size_t& order_u)
{
    std::vector<std::size_t> orders(1);
    orders[0] = order_u;
    return rDummy.CreateRegularPatch<1>(orders);
}

NURBSPatch<2>::Pointer NURBSPatchLibrary_CreateRectangularPatch(NURBSPatchLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v)
{
    std::vector<std::size_t> orders(2);
    orders[0] = order_u;
    orders[1] = order_v;
    return rDummy.CreateRegularPatch<2>(orders);
}

NURBSPatch<3>::Pointer NURBSPatchLibrary_CreateCubicPatch(NURBSPatchLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v, const std::size_t& order_w)
{
    std::vector<std::size_t> orders(3);
    orders[0] = order_u;
    orders[1] = order_v;
    orders[2] = order_w;
    return rDummy.CreateRegularPatch<3>(orders);
}

void NURBSPatchUtility_ExportGeo1(NURBSPatchUtility& rDummy, NURBSPatch<1>::Pointer pPatch, const std::string& filename)
{
    rDummy.ExportGeo(pPatch, filename);
}

void NURBSPatchUtility_ExportGeo2(NURBSPatchUtility& rDummy, NURBSPatch<2>::Pointer pPatch, const std::string& filename)
{
    rDummy.ExportGeo(pPatch, filename);
}

void NURBSPatchUtility_ExportGeo3(NURBSPatchUtility& rDummy, NURBSPatch<3>::Pointer pPatch, const std::string& filename)
{
    rDummy.ExportGeo(pPatch, filename);
}

GridFunction<1, ControlPoint<double> >::Pointer GridFunctionUtility_CreateLinearControlPointGrid(
        GridFunctionUtility& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u,
        const double& end_x, const double& end_y, const double& end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateRegularControlPointGrid<1>(start, ngrid, end);
}

GridFunction<2, ControlPoint<double> >::Pointer GridFunctionUtility_CreateRectangularControlPointGrid(
        GridFunctionUtility& rDummy,
        const double& start_x, const double& start_y,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& end_x, const double& end_y)
{
    std::vector<double> start(2);
    start[0] = start_x;
    start[1] = start_y;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> end(2);
    end[0] = end_x;
    end[1] = end_y;

    return rDummy.CreateRegularControlPointGrid<2>(start, ngrid, end);
}

GridFunction<3, ControlPoint<double> >::Pointer GridFunctionUtility_CreateCubicControlPointGrid(
        GridFunctionUtility& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        const double& end_x, const double& end_y, const double& end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateRegularControlPointGrid<3>(start, ngrid, end);
}

template<std::size_t TDim>
void MultiPatchRefinementUtility_InsertKnots(MultiPatchRefinementUtility& rDummy,
        typename NURBSPatch<TDim>::Pointer& pPatch,
        boost::python::list ins_knots)
{
    std::vector<std::vector<double> > ins_knots_array(TDim);
    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& ins_knots_x,
                std::make_pair(iterator_value_type(ins_knots), // begin
                iterator_value_type() ) ) // end
    {
        std::vector<double> knots;

        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        BOOST_FOREACH(const iterator_value_type2::value_type& knot,
                    std::make_pair(iterator_value_type2(ins_knots_x), // begin
                    iterator_value_type2() ) ) // end
        {
            knots.push_back(knot);
        }

        ins_knots_array[dim++] = knots;
    }

    rDummy.InsertKnots<TDim>(pPatch, ins_knots_array);
}

template<std::size_t TDim>
void MultiPatchRefinementUtility_DegreeElevate(MultiPatchRefinementUtility& rDummy,
        typename NURBSPatch<TDim>::Pointer& pPatch,
        boost::python::list order_increment)
{
    std::vector<std::size_t> order_incr_array(TDim);
    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& t,
                std::make_pair(iterator_value_type(order_increment), // begin
                iterator_value_type() ) ) // end
    {
        order_incr_array[dim++] = static_cast<std::size_t>(t);
    }

    rDummy.DegreeElevate<TDim>(pPatch, order_incr_array);
}

template<std::size_t TDim>
void MultiPatch_MakeNeighbor(MultiPatch<TDim>& rDummy, typename Patch<TDim>::Pointer pPatch1, BoundarySide side1,
            typename Patch<TDim>::Pointer pPatch2, BoundarySide side2)
{
    rDummy.MakeNeighbor(pPatch1, side1, pPatch2, side2);
}

/////////////TESTING//////////////////
template<std::size_t TDim>
void MultiPatch_Test_MakeNeighbor(MultiPatch<TDim>& rDummy)
{
    std::vector<std::size_t> orders(TDim);

    if (TDim == 2)
    {
        orders[0] = 3;
        orders[1] = 3;
    }
    else if (TDim == 3)
    {
        orders[0] = 3;
        orders[1] = 3;
        orders[2] = 3;
    }

    typename Patch<TDim>::Pointer pPatch1 = NURBSPatchLibrary::CreateRegularPatch<TDim>(orders);
    typename Patch<TDim>::Pointer pPatch2 = NURBSPatchLibrary::CreateRegularPatch<TDim>(orders);

    rDummy.AddPatch(pPatch1);
    rDummy.AddPatch(pPatch2);
    rDummy.ResetId();

    // rDummy.MakeNeighbor(pPatch1, _RIGHT_, pPatch2, _LEFT_);
    MultiPatch_MakeNeighbor<TDim>(rDummy, pPatch1, _RIGHT_, pPatch2, _LEFT_);

    // KRATOS_WATCH(*pPatch1)
    // KRATOS_WATCH(*pPatch2)
}
//////////////////////////////////////

template<class TPatchType>
TPatchType& GetReference(typename TPatchType::Pointer& dummy)
{
    return *dummy;
}

template<class TPatchType>
typename TPatchType::Pointer Patch_pGetLeft(typename TPatchType::Pointer dummy)
{
    return dummy->pLeft();
}

template<class TPatchType>
typename TPatchType::Pointer Patch_pGetRight(typename TPatchType::Pointer dummy)
{
    return dummy->pRight();
}

template<std::size_t TDim>
void IsogeometricApplication_AddPatchesToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim> >
    // class_<Patch<TDim>, typename Patch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&>())
    .def("SetControlPointGrid", &Patch<TDim>::SetControlPointGrid)
    .def("Order", &Patch<TDim>::Order)
    .def("TotalNumber", &Patch<TDim>::TotalNumber)
    .def("Left", &Patch_pGetLeft<Patch<TDim> >)
    .def("Right", &Patch_pGetRight<Patch<TDim> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Patch" << TDim << "DPointer";
    class_<typename Patch<TDim>::Pointer>
    (ss.str().c_str(), init<typename Patch<TDim>::Pointer>())
    .def("GetReference", GetReference<Patch<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "NURBSPatch" << TDim << "D";
    class_<NURBSPatch<TDim>, bases<Patch<TDim> > >
    // class_<NURBSPatch<TDim>, typename NURBSPatch<TDim>::Pointer, bases<Patch<TDim> > >
    (ss.str().c_str(), init<const std::size_t&>())
    .def("Number", &NURBSPatch<TDim>::Number)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "NURBSPatch" << TDim << "DPointer";
    class_<typename NURBSPatch<TDim>::Pointer>
    (ss.str().c_str(), init<typename NURBSPatch<TDim>::Pointer>())
    .def("GetReference", GetReference<NURBSPatch<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "HnMesh" << TDim << "D";
    class_<HnMesh<TDim>, bases<Patch<TDim> > >
    // class_<HnMesh<TDim>, typename HnMesh<TDim>::Pointer, bases<Patch<TDim> > >
    (ss.str().c_str(), init<const std::size_t&, const std::string&>())
    .def("SetEchoLevel", &HnMesh<TDim>::SetEchoLevel)
    .def("ReadMesh", &HnMesh<TDim>::ReadMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("Refine", &HnMesh<TDim>::Refine) // use this for debugging only, use RefineNodes and LinearDependencyRefine instead
    .def("RefineNodes", &HnMesh<TDim>::RefineNodes)
    .def("LinearDependencyRefine", &HnMesh<TDim>::LinearDependencyRefine)
    .def("BuildMesh", &HnMesh<TDim>::BuildMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportCellTopology", &HnMesh<TDim>::ExportCellTopology)
    .def("ExportCellGeology", &HnMesh<TDim>::ExportCellGeology)
//    .def("ExportRefinedDomain", &HnMesh<TDim>::ExportRefinedDomain)
    .def("ExportSupportDomain", &HnMesh<TDim>::ExportSupportDomain)
    .def("ExportMatlab", &HnMesh<TDim>::ExportMatlab)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportMDPA", &HnMesh<TDim>::ExportMDPA)
    .def("ExportMDPA2", &HnMesh<TDim>::ExportMDPA2)
    .def("ExportPostMDPA", &HnMesh<TDim>::ExportPostMDPA)
    .def("ExportCellGeologyAsPostMDPA", &HnMesh<TDim>::ExportCellGeologyAsPostMDPA)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("PrintKnotVectors", &HnMesh<TDim>::PrintKnotVectors)
    .def("PrintCells", &HnMesh<TDim>::PrintCells)
    .def("PrintBasisFuncs", &HnMesh<TDim>::PrintBasisFuncs)
    .def("PrintRefinementHistory", &HnMesh<TDim>::PrintRefinementHistory)
    .def("CheckNestedSpace", &HnMesh<TDim>::CheckNestedSpace)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "HnMesh" << TDim << "DPointer";
    class_<typename HnMesh<TDim>::Pointer>
    (ss.str().c_str(), init<typename HnMesh<TDim>::Pointer>())
    .def("GetReference", GetReference<HnMesh<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiPatch" << TDim << "D";
    class_<MultiPatch<TDim>, typename MultiPatch<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("ResetId", &MultiPatch<TDim>::ResetId)
    .def("AddPatch", &MultiPatch<TDim>::AddPatch)
    .def("GetPatch", &MultiPatch<TDim>::GetPatch)
    .def("MakeNeighbor", &MultiPatch_MakeNeighbor<TDim>)
    .def("Test_MakeNeighbor", &MultiPatch_Test_MakeNeighbor<TDim>)
    .def("PrintAddress", &MultiPatch<TDim>::PrintAddress)
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddGridFunctions()
{
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<GridFunction<0, ControlPoint<double> >, GridFunction<0, ControlPoint<double> >::Pointer, boost::noncopyable>
    ("ControlPointGridFunction", init<>())
    .def("Size", &GridFunction<0, ControlPoint<double> >::Size)
    .def(self_ns::str(self))
    ;

    class_<GridFunction<0, double>, GridFunction<0, double>::Pointer, boost::noncopyable>
    ("DoubleGridFunction", init<>())
    .def("Size", &GridFunction<0, double>::Size)
    .def(self_ns::str(self))
    ;

    class_<GridFunction<0, array_1d<double, 3> >, GridFunction<0, array_1d<double, 3> >::Pointer, boost::noncopyable>
    ("Array1DGridFunction", init<>())
    .def("Size", &GridFunction<0, array_1d<double, 3> >::Size)
    .def(self_ns::str(self))
    ;

    class_<GridFunction<0, Vector>, GridFunction<0, Vector>::Pointer, boost::noncopyable>
    ("VectorGridFunction", init<>())
    .def("Size", &GridFunction<0, Vector>::Size)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<GridFunction<1, ControlPoint<double> >, GridFunction<1, ControlPoint<double> >::Pointer, bases<GridFunction<0, ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointGridFunction1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<1, double>, GridFunction<1, double>::Pointer, bases<GridFunction<0, double> >, boost::noncopyable>
    ("DoubleGridFunction1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<1, array_1d<double, 3> >, GridFunction<1, array_1d<double, 3> >::Pointer, bases<GridFunction<0, array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DGridFunction1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<1, Vector>, GridFunction<1, Vector>::Pointer, bases<GridFunction<0, Vector> >, boost::noncopyable>
    ("VectorGridFunction1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<GridFunction<2, ControlPoint<double> >, GridFunction<2, ControlPoint<double> >::Pointer, bases<GridFunction<0, ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointGridFunction2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<2, double>, GridFunction<2, double>::Pointer, bases<GridFunction<0, double> >, boost::noncopyable>
    ("DoubleGridFunction2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<2, array_1d<double, 3> >, GridFunction<2, array_1d<double, 3> >::Pointer, bases<GridFunction<0, array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DGridFunction2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<2, Vector>, GridFunction<2, Vector>::Pointer, bases<GridFunction<0, Vector> >, boost::noncopyable>
    ("VectorGridFunction2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<GridFunction<3, ControlPoint<double> >, GridFunction<3, ControlPoint<double> >::Pointer, bases<GridFunction<0, ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointGridFunction3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<3, double>, GridFunction<3, double>::Pointer, bases<GridFunction<0, double> >, boost::noncopyable>
    ("DoubleGridFunction3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<3, array_1d<double, 3> >, GridFunction<3, array_1d<double, 3> >::Pointer, bases<GridFunction<0, array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DGridFunction3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<GridFunction<3, Vector>, GridFunction<3, Vector>::Pointer, bases<GridFunction<0, Vector> >, boost::noncopyable>
    ("VectorGridFunction3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddCustomUtilities2ToPython()
{
    /////////////////////////////////////////////////////////////////
    ///////////////////////NURBS/////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<DomainManager, DomainManager::Pointer, boost::noncopyable>
    ("DomainManager", init<std::size_t>())
    ;

    class_<DomainManager2D, DomainManager2D::Pointer, boost::noncopyable>
    ("DomainManager2D", init<std::size_t>())
    .def("AddXcoord", &DomainManager2D::AddXcoord)
    .def("AddYcoord", &DomainManager2D::AddYcoord)
    .def("AddCell", &DomainManager2D::AddCell)
    .def("IsInside", &DomainManager2D::IsInside)
    .def(self_ns::str(self))
    ;

    class_<NURBSPatchLibrary, NURBSPatchLibrary::Pointer, boost::noncopyable>
    ("NURBSPatchLibrary", init<>())
    .def("CreateLinearPatch", &NURBSPatchLibrary_CreateLinearPatch)
    .def("CreateRectangularPatch", &NURBSPatchLibrary_CreateRectangularPatch)
    .def("CreateCubicPatch", &NURBSPatchLibrary_CreateCubicPatch)
    ;

    class_<NURBSPatchUtility, NURBSPatchUtility::Pointer, boost::noncopyable>
    ("NURBSPatchUtility", init<>())
    .def("ExportGeo", &NURBSPatchUtility_ExportGeo1)
    .def("ExportGeo", &NURBSPatchUtility_ExportGeo2)
    .def("ExportGeo", &NURBSPatchUtility_ExportGeo3)
    ;

    class_<GridFunctionUtility, GridFunctionUtility::Pointer, boost::noncopyable>
    ("GridFunctionUtility", init<>())
    .def("CreateLinearControlPointGrid", &GridFunctionUtility_CreateLinearControlPointGrid)
    .def("CreateRectangularControlPointGrid", &GridFunctionUtility_CreateRectangularControlPointGrid)
    .def("CreateCubicControlPointGrid", &GridFunctionUtility_CreateCubicControlPointGrid)
    ;

    class_<MultiPatchRefinementUtility, MultiPatchRefinementUtility::Pointer, boost::noncopyable>
    ("MultiPatchRefinementUtility", init<>())
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<1>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<2>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<3>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<1>)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////GRIDS/////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddGridFunctions();

    /////////////////////////////////////////////////////////////////
    ///////////////////////TSPLINES//////////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<TsMesh2D, TsMesh2D::Pointer, boost::noncopyable>
    ("TsMesh2D", init<>())
    .def("BeginConstruct", &TsMesh2D::BeginConstruct)
    .def("EndConstruct", &TsMesh2D::EndConstruct)
    .def("ReadFromFile", &TsMesh2D::ReadFromFile)
    .def("ExportMatlab", &TsMesh2D::ExportMatlab)
    .def("BuildExtendedTmesh", &TsMesh2D::BuildExtendedTmesh)
    .def("IsAnalysisSuitable", &TsMesh2D::IsAnalysisSuitable)
    .def("BuildAnchors", &TsMesh2D::BuildAnchors)
    .def("BuildCells", &TsMesh2D::BuildCells)
    .def("ExportMDPA", &TsMesh2D::ExportMDPA)
//    .def("FindKnots2", &TsMesh2D::FindKnots2)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////HIERARCHICAL NURBS////////////////////////
    /////////////////////////////////////////////////////////////////

    enum_<HN_ECHO_FLAGS>("HN_ECHO_FLAGS")
    .value("ECHO_REFIMENT", ECHO_REFIMENT)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////PATCHES///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    enum_<BoundarySide>("BoundarySide")
    .value("Left", _LEFT_)
    .value("Right", _RIGHT_)
    .value("Top", _TOP_)
    .value("Bottom", _BOTTOM_)
    .value("Front", _FRONT_)
    .value("Back", _BACK_)
    ;

    IsogeometricApplication_AddPatchesToPython<1>();
    IsogeometricApplication_AddPatchesToPython<2>();
    IsogeometricApplication_AddPatchesToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos


/*
LICENSE: see isogeometric_application/LICENSE.txt
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
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/trans/translation.h"
#include "custom_utilities/trans/rotation.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/nurbs/regular_control_grid.h"
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/multipatch_refinement_utility.h"
#include "custom_utilities/hierarchical_bsplines/hb_mesh.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"
#include "custom_python/add_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

double ControlPoint_GetWX(ControlPoint<double>& rDummy)
{
    return rDummy.WX();
}

void ControlPoint_SetWX(ControlPoint<double>& rDummy, const double& newWX)
{
    rDummy.WX() = newWX;
}

double ControlPoint_GetWY(ControlPoint<double>& rDummy)
{
    return rDummy.WY();
}

void ControlPoint_SetWY(ControlPoint<double>& rDummy, const double& newWY)
{
    rDummy.WY() = newWY;
}

double ControlPoint_GetWZ(ControlPoint<double>& rDummy)
{
    return rDummy.WZ();
}

void ControlPoint_SetWZ(ControlPoint<double>& rDummy, const double& newWZ)
{
    rDummy.WZ() = newWZ;
}

double ControlPoint_GetW(ControlPoint<double>& rDummy)
{
    return rDummy.W();
}

void ControlPoint_SetW(ControlPoint<double>& rDummy, const double& newW)
{
    rDummy.W() = newW;
}

void ControlPoint_ApplyTransformation(ControlPoint<double>& rDummy, const Transformation<double>& trans)
{
    rDummy.ApplyTransformation(trans);
}

////////////////////////////////////////

template<typename TDataType>
TDataType ControlGrid_GetItem(ControlGrid<TDataType>& rDummy, int index)
{
    return rDummy.GetData(index);
}

template<typename TDataType>
void ControlGrid_SetItem(ControlGrid<TDataType>& rDummy, int index, const TDataType& value)
{
    rDummy.SetData(index, value);
}

////////////////////////////////////////

ControlGrid<ControlPoint<double> >::Pointer ControlGridUtility_CreateLinearControlPointGrid(
        ControlGridUtility& rDummy,
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

ControlGrid<ControlPoint<double> >::Pointer ControlGridUtility_CreateRectangularControlPointGrid1(
        ControlGridUtility& rDummy,
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

ControlGrid<ControlPoint<double> >::Pointer ControlGridUtility_CreateRectangularControlPointGrid2(
        ControlGridUtility& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& space1_x, const double& space1_y, const double& space1_z,
        const double& space2_x, const double& space2_y, const double& space2_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> space1(3);
    space1[0] = space1_x;
    space1[1] = space1_y;
    space1[2] = space1_z;

    std::vector<double> space2(3);
    space2[0] = space2_x;
    space2[1] = space2_y;
    space2[2] = space2_z;

    std::vector<std::vector<double> > spacing_vectors(2);
    spacing_vectors[0] = space1;
    spacing_vectors[1] = space2;

    return rDummy.CreateRegularControlPointGrid<2>(start, ngrid, spacing_vectors);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridUtility_CreateCubicControlPointGrid1(
        ControlGridUtility& rDummy,
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

ControlGrid<ControlPoint<double> >::Pointer ControlGridUtility_CreateCubicControlPointGrid2(
        ControlGridUtility& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        boost::python::list spacing_vectors_data)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<std::vector<double> > spacing_vectors;
    std::size_t cnt1 = 0, cnt2 = 0;
    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& vect, std::make_pair(iterator_value_type(spacing_vectors_data), iterator_value_type() ) )
    {
        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        std::vector<double> space_vect;
        BOOST_FOREACH(const iterator_value_type2::value_type& v, std::make_pair(iterator_value_type2(vect), iterator_value_type2() ) )
        {
            space_vect.push_back(v);
        }
        spacing_vectors.push_back(space_vect);
    }

    return rDummy.CreateRegularControlPointGrid<3>(start, ngrid, spacing_vectors);
}

////////////////////////////////////////

BSplinesFESpace<1>::Pointer BSplinesFESpaceLibrary_CreateLinearFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u)
{
    std::vector<std::size_t> orders(1);
    orders[0] = order_u;
    return rDummy.CreateRegularFESpace<1>(orders);
}

BSplinesFESpace<2>::Pointer BSplinesFESpaceLibrary_CreateRectangularFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v)
{
    std::vector<std::size_t> orders(2);
    orders[0] = order_u;
    orders[1] = order_v;
    return rDummy.CreateRegularFESpace<2>(orders);
}

BSplinesFESpace<3>::Pointer BSplinesFESpaceLibrary_CreateCubicFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v, const std::size_t& order_w)
{
    std::vector<std::size_t> orders(3);
    orders[0] = order_u;
    orders[1] = order_v;
    orders[2] = order_w;
    return rDummy.CreateRegularFESpace<3>(orders);
}

////////////////////////////////////////

template<int TDim, typename TDataType>
typename FESpace<TDim>::Pointer GridFunction_GetFESpace(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pFESpace();
}

template<int TDim, typename TDataType>
void GridFunction_SetFESpace(GridFunction<TDim, TDataType>& rDummy, typename FESpace<TDim>::Pointer pNewFESpace)
{
    rDummy.SetFESpace(pNewFESpace);
}

template<int TDim, typename TDataType>
typename ControlGrid<TDataType>::Pointer GridFunction_GetControlGrid(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pControlGrid();
}

template<int TDim, typename TDataType>
void GridFunction_SetControlGrid(GridFunction<TDim, TDataType>& rDummy, typename ControlGrid<TDataType>::Pointer pNewControlGrid)
{
    rDummy.SetControlGrid(pNewControlGrid);
}

////////////////////////////////////////

template<class TPatchType>
TPatchType& GetReference(typename TPatchType::Pointer& dummy)
{
    return *dummy;
}

template<class TPatchType>
std::size_t Patch_GetId(TPatchType& rDummy)
{
    return rDummy.Id();
}

template<class TPatchType>
void Patch_SetId(TPatchType& rDummy, const std::size_t& Id)
{
    rDummy.SetId(Id);
}

template<class TPatchType>
typename TPatchType::Pointer Patch_pGetNeighbor(TPatchType& rDummy, const BoundarySide& side)
{
    return rDummy.pNeighbor(side);
}

template<class TPatchType, class TMultiPatchType>
typename TPatchType::Pointer MultiPatch_GetItem(TMultiPatchType& rDummy, std::size_t index)
{
    return rDummy.pGetPatch(index);
}

template<int TDim>
void MultiPatch_MakeNeighbor(MultiPatch<TDim>& rDummy, typename Patch<TDim>::Pointer pPatch1, BoundarySide side1,
           typename Patch<TDim>::Pointer pPatch2, BoundarySide side2)
{
   rDummy.MakeNeighbor(pPatch1, side1, pPatch2, side2);
}

template<int TDim>
std::size_t MultiPatch_Enumerate(MultiPatch<TDim>& rDummy)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate();
    KRATOS_WATCH(system_size)

    return system_size;
}

template<int TDim>
void MultiPatchRefinementUtility_InsertKnots(MultiPatchRefinementUtility& rDummy,
       typename Patch<TDim>::Pointer& pPatch,
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

template<int TDim>
void MultiPatchRefinementUtility_DegreeElevate(MultiPatchRefinementUtility& rDummy,
       typename Patch<TDim>::Pointer& pPatch,
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

template<int TDim, class TExporter, class TPatchType>
void MultiPatchExporter_Export(TExporter& rDummy,
        typename TPatchType::Pointer pPatch, const std::string& filename)
{
    rDummy.template Export<TDim>(pPatch, filename);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void IsogeometricApplication_AddTransformation()
{
    class_<Transformation<double>, Transformation<double>::Pointer, boost::noncopyable>
    ("Transformation", init<>())
    .def("AppendTransformation", &Transformation<double>::AppendTransformation)
    // .def(boost::python::operators<boost::python::op_mul>());
    .def(self_ns::str(self))
    ;

    class_<Translation<double>, Translation<double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("Translation", init<const double&, const double&, const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<0, double>, Rotation<0, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationX", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<1, double>, Rotation<1, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationY", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<2, double>, Rotation<2, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationZ", init<const double&>())
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddControlPoint()
{
    class_<ControlPoint<double>, ControlPoint<double>::Pointer>
    ("ControlPoint", init<>())
    .add_property("WX", ControlPoint_GetWX, ControlPoint_SetWX)
    .add_property("WY", ControlPoint_GetWY, ControlPoint_SetWY)
    .add_property("WZ", ControlPoint_GetWZ, ControlPoint_SetWZ)
    .add_property("W", ControlPoint_GetW, ControlPoint_SetW)
    .def("ApplyTransformation", &ControlPoint_ApplyTransformation)
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddControlGrids()
{
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<ControlGrid<ControlPoint<double> >, ControlGrid<ControlPoint<double> >::Pointer, boost::noncopyable>
    ("ControlPointControlGrid", init<>())
    .def("Size", &ControlGrid<ControlPoint<double> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<ControlPoint<double> >)
    .def("__getitem__", &ControlGrid_GetItem<ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<double>, ControlGrid<double>::Pointer, boost::noncopyable>
    ("DoubleControlGrid", init<>())
    .def("Size", &ControlGrid<double>::Size)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<array_1d<double, 3> >, ControlGrid<array_1d<double, 3> >::Pointer, boost::noncopyable>
    ("Array1DControlGrid", init<>())
    .def("Size", &ControlGrid<array_1d<double, 3> >::Size)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<Vector>, ControlGrid<Vector>::Pointer, boost::noncopyable>
    ("VectorControlGrid", init<>())
    .def("Size", &ControlGrid<Vector>::Size)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<RegularControlGrid<1, ControlPoint<double> >, RegularControlGrid<1, ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointControlGrid1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<1, double>, RegularControlGrid<1, double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("DoubleControlGrid1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<1, array_1d<double, 3> >, RegularControlGrid<1, array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DControlGrid1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<1, Vector>, RegularControlGrid<1, Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("VectorControlGrid1D", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<RegularControlGrid<2, ControlPoint<double> >, RegularControlGrid<2, ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<2, double>, RegularControlGrid<2, double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("DoubleControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<2, array_1d<double, 3> >, RegularControlGrid<2, array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<2, Vector>, RegularControlGrid<2, Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("VectorControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<RegularControlGrid<3, ControlPoint<double> >, RegularControlGrid<3, ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("ControlPointControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<3, double>, RegularControlGrid<3, double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("DoubleControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<3, array_1d<double, 3> >, RegularControlGrid<3, array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("Array1DControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<RegularControlGrid<3, Vector>, RegularControlGrid<3, Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("VectorControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddFESpacesToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "FESpace" << TDim << "D";
    class_<FESpace<TDim>, typename FESpace<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Order", &FESpace<TDim>::Order)
    .def("TotalNumber", &FESpace<TDim>::TotalNumber)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "BSplinesFESpace" << TDim << "D";
    class_<BSplinesFESpace<TDim>, typename BSplinesFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Number", &BSplinesFESpace<TDim>::Number)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddGridFunctionsToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "ControlPointGridFunction" << TDim << "D";
    class_<GridFunction<TDim, ControlPoint<double> >, typename GridFunction<TDim, ControlPoint<double> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<ControlPoint<double> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, ControlPoint<double> >, GridFunction_SetFESpace<TDim, ControlPoint<double> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, ControlPoint<double> >, GridFunction_SetControlGrid<TDim, ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "DoubleGridFunction" << TDim << "D";
    class_<GridFunction<TDim, double>, typename GridFunction<TDim, double>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<double>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, double>, GridFunction_SetFESpace<TDim, double>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, double>, GridFunction_SetControlGrid<TDim, double>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Array1DGridFunction" << TDim << "D";
    class_<GridFunction<TDim, array_1d<double, 3> >, typename GridFunction<TDim, array_1d<double, 3> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<array_1d<double, 3> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, array_1d<double, 3> >, GridFunction_SetFESpace<TDim, array_1d<double, 3> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, array_1d<double, 3> >, GridFunction_SetControlGrid<TDim, array_1d<double, 3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "VectorGridFunction" << TDim << "D";
    class_<GridFunction<TDim, Vector>, typename GridFunction<TDim, Vector>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<Vector>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, Vector>, GridFunction_SetFESpace<TDim, Vector>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, Vector>, GridFunction_SetControlGrid<TDim, Vector>)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddPatchesToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim> >
    // class_<Patch<TDim>, typename Patch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&, typename FESpace<TDim>::Pointer>())
    .add_property("Id", &Patch_GetId<Patch<TDim> >, &Patch_SetId<Patch<TDim> >)
    .def("CreateControlPointGridFunction", &Patch<TDim>::CreateControlPointGridFunction)
    .def("CreateDoubleGridFunction", &Patch<TDim>::CreateDoubleGridFunction)
    .def("CreateArray1DGridFunction", &Patch<TDim>::CreateArray1DGridFunction)
    .def("CreateVectorGridFunction", &Patch<TDim>::CreateVectorGridFunction)
    .def("ApplyTransformation", &Patch<TDim>::ApplyTransformation)
    .def("Order", &Patch<TDim>::Order)
    .def("TotalNumber", &Patch<TDim>::TotalNumber)
    .def("Neighbor", &Patch_pGetNeighbor<Patch<TDim> >)
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
    ss << "HBMesh" << TDim << "D";
    class_<HBMesh<TDim>, bases<Patch<TDim> > >
    // class_<HBMesh<TDim>, typename HBMesh<TDim>::Pointer, bases<Patch<TDim> > >
    (ss.str().c_str(), init<const std::size_t&, const std::string&>())
    .def("SetEchoLevel", &HBMesh<TDim>::SetEchoLevel)
    .def("ReadMesh", &HBMesh<TDim>::ReadMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("Refine", &HBMesh<TDim>::Refine) // use this for debugging only, use RefineNodes and LinearDependencyRefine instead
    .def("RefineNodes", &HBMesh<TDim>::RefineNodes)
    .def("LinearDependencyRefine", &HBMesh<TDim>::LinearDependencyRefine)
    .def("BuildMesh", &HBMesh<TDim>::BuildMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportCellTopology", &HBMesh<TDim>::ExportCellTopology)
    .def("ExportCellGeology", &HBMesh<TDim>::ExportCellGeology)
    //    .def("ExportRefinedDomain", &HBMesh<TDim>::ExportRefinedDomain)
    .def("ExportSupportDomain", &HBMesh<TDim>::ExportSupportDomain)
    .def("ExportMatlab", &HBMesh<TDim>::ExportMatlab)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportMDPA", &HBMesh<TDim>::ExportMDPA)
    .def("ExportMDPA2", &HBMesh<TDim>::ExportMDPA2)
    .def("ExportPostMDPA", &HBMesh<TDim>::ExportPostMDPA)
    .def("ExportCellGeologyAsPostMDPA", &HBMesh<TDim>::ExportCellGeologyAsPostMDPA)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("PrintKnotVectors", &HBMesh<TDim>::PrintKnotVectors)
    .def("PrintCells", &HBMesh<TDim>::PrintCells)
    .def("PrintBasisFuncs", &HBMesh<TDim>::PrintBasisFuncs)
    .def("PrintRefinementHistory", &HBMesh<TDim>::PrintRefinementHistory)
    .def("CheckNestedSpace", &HBMesh<TDim>::CheckNestedSpace)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "HBMesh" << TDim << "DPointer";
    class_<typename HBMesh<TDim>::Pointer>
    (ss.str().c_str(), init<typename HBMesh<TDim>::Pointer>())
    .def("GetReference", GetReference<HBMesh<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiPatch" << TDim << "D";
    class_<MultiPatch<TDim>, typename MultiPatch<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("ResetId", &MultiPatch<TDim>::ResetId)
    .def("AddPatch", &MultiPatch<TDim>::AddPatch)
    .def("__getitem__", &MultiPatch_GetItem<Patch<TDim>, MultiPatch<TDim> >)
    .def("MakeNeighbor", &MultiPatch_MakeNeighbor<TDim>)
    .def("EquationSystemSize", &MultiPatch<TDim>::EquationSystemSize)
    .def("Enumerate", &MultiPatch_Enumerate<TDim>)
    .def("PrintAddress", &MultiPatch<TDim>::PrintAddress)
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddImportExportToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "MultiNURBSPatchGeoExporter";
    class_<MultiNURBSPatchGeoExporter, MultiNURBSPatchGeoExporter::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiNURBSPatchMatlabExporter";
    class_<MultiNURBSPatchMatlabExporter, MultiNURBSPatchMatlabExporter::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiNURBSPatchGLVisExporter";
    class_<MultiNURBSPatchGLVisExporter, MultiNURBSPatchGLVisExporter::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGLVisExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGLVisExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

}

void IsogeometricApplication_AddCustomUtilities2ToPython()
{
    /////////////////////////////////////////////////////////////////
    ///////////////////////SUPPORT DOMAIN////////////////////////////
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

    /////////////////////////////////////////////////////////////////
    ///////////////////TRANSFORMATION////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddTransformation();

    /////////////////////////////////////////////////////////////////
    ///////////////////CONTROL POINT/////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlPoint();

    /////////////////////////////////////////////////////////////////
    ///////////////////////CONTROL GRIDS/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlGrids();

    class_<ControlGridUtility, ControlGridUtility::Pointer, boost::noncopyable>
    ("ControlGridUtility", init<>())
    .def("CreateLinearControlPointGrid", &ControlGridUtility_CreateLinearControlPointGrid)
    .def("CreateRectangularControlPointGrid", &ControlGridUtility_CreateRectangularControlPointGrid1)
    .def("CreateRectangularControlPointGrid", &ControlGridUtility_CreateRectangularControlPointGrid2)
    .def("CreateCubicControlPointGrid", &ControlGridUtility_CreateCubicControlPointGrid1)
    .def("CreateCubicControlPointGrid", &ControlGridUtility_CreateCubicControlPointGrid2)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////FESpace///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddFESpacesToPython<1>();
    IsogeometricApplication_AddFESpacesToPython<2>();
    IsogeometricApplication_AddFESpacesToPython<3>();

    class_<BSplinesFESpaceLibrary, BSplinesFESpaceLibrary::Pointer, boost::noncopyable>
    ("BSplinesFESpaceLibrary", init<>())
    .def("CreateLinearFESpace", &BSplinesFESpaceLibrary_CreateLinearFESpace)
    .def("CreateRectangularFESpace", &BSplinesFESpaceLibrary_CreateRectangularFESpace)
    .def("CreateCubicFESpace", &BSplinesFESpaceLibrary_CreateCubicFESpace)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////Grid Functions////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddGridFunctionsToPython<1>();
    IsogeometricApplication_AddGridFunctionsToPython<2>();
    IsogeometricApplication_AddGridFunctionsToPython<3>();

    /////////////////////////////////////////////////////////////////
    ///////////////////////GRID FUNCTIONS////////////////////////////
    /////////////////////////////////////////////////////////////////


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
    ///////////////////////HIERARCHICAL BSplines/////////////////////
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

    class_<MultiPatchUtility, MultiPatchUtility::Pointer, boost::noncopyable>
    ("MultiPatchUtility", init<>())
    .def("CreatePatchPointer", &MultiPatchUtility::CreatePatchPointer<1>)
    .def("CreatePatchPointer", &MultiPatchUtility::CreatePatchPointer<2>)
    .def("CreatePatchPointer", &MultiPatchUtility::CreatePatchPointer<3>)
    ;

    class_<MultiPatchRefinementUtility, MultiPatchRefinementUtility::Pointer, boost::noncopyable>
    ("MultiPatchRefinementUtility", init<>())
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<1>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<2>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<3>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<1>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<2>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<3>)
    ;

    class_<BSplinesPatchUtility, BSplinesPatchUtility::Pointer, boost::noncopyable>
    ("BSplinesPatchUtility", init<>())
    .def("CreateConnectedPatch", &BSplinesPatchUtility::CreateConnectedPatch<2>)
    .def("CreateConnectedPatch", &BSplinesPatchUtility::CreateConnectedPatch<3>)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////IMPORT/EXPORT/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddImportExportToPython();

}

}  // namespace Python.

} // Namespace Kratos


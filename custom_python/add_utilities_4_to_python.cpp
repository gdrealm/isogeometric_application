/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 21, 2017 $
//   Revision:            $Revision: 1.0 $
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
#include "custom_utilities/nonconforming_multipatch_lagrange_mesh.h"
#include "custom_utilities/multipatch_model_part.h"
#include "custom_utilities/multi_multipatch_model_part.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<class T>
ModelPart& MultiPatchModelPart_GetModelPart(T& rDummy)
{
    return *(rDummy.pModelPart());
}

template<class T>
typename T::MultiPatchType& MultiPatchModelPart_GetMultiPatch(T& rDummy)
{
    return *(rDummy.pMultiPatch());
}

template<class T>
typename T::MultiPatchType& MultiPatchModelPart_GetMultiPatch2(T& rDummy, const std::size_t& i)
{
    return *(rDummy.pMultiPatch(i));
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddMeshToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "NonConformingMultipatchLagrangeMesh" << TDim << "D";
    class_<NonConformingMultipatchLagrangeMesh<TDim>, typename NonConformingMultipatchLagrangeMesh<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("SetBaseElementName", &NonConformingMultipatchLagrangeMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastElemId)
    .def("SetLastPropId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastPropId)
    .def("SetDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetDivision)
    .def("SetUniformDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetUniformDivision)
    .def("WriteModelPart", &NonConformingMultipatchLagrangeMesh<TDim>::WriteModelPart)
    .def(self_ns::str(self))
    ;

}

template<int TDim>
void IsogeometricApplication_AddModelPartToPython()
{
    std::stringstream ss;

    typedef MultiPatchModelPart<TDim> MultiPatchModelPartType;
    ss.str(std::string());
    ss << "MultiPatchModelPart" << TDim << "D";
    class_<MultiPatchModelPartType, typename MultiPatchModelPartType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("BeginModelPart", &MultiPatchModelPartType::BeginModelPart)
    .def("CreateNodes", &MultiPatchModelPartType::CreateNodes)
    .def("AddElements", &MultiPatchModelPartType::AddElements)
    .def("AddConditions", &MultiPatchModelPartType::AddConditions)
    .def("EndModelPart", &MultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiPatchModelPartType>, return_internal_reference<>())
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch<MultiPatchModelPartType>, return_internal_reference<>())
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    .def(self_ns::str(self))
    ;

    typedef MultiMultiPatchModelPart<TDim> MultiMultiPatchModelPartType;
    ss.str(std::string());
    ss << "MultiMultiPatchModelPart" << TDim << "D";
    class_<MultiMultiPatchModelPartType, typename MultiMultiPatchModelPartType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("AddMultiPatch", &MultiMultiPatchModelPartType::AddMultiPatch)
    .def("BeginModelPart", &MultiMultiPatchModelPartType::BeginModelPart)
    .def("CreateNodes", &MultiMultiPatchModelPartType::CreateNodes)
    // .def("AddElements", &MultiMultiPatchModelPartType::AddElements)
    // .def("AddConditions", &MultiMultiPatchModelPartType::AddConditions)
    .def("EndModelPart", &MultiMultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiMultiPatchModelPartType>, return_internal_reference<>())
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch2<MultiMultiPatchModelPartType>, return_internal_reference<>())
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    .def(self_ns::str(self))
    ;
}


void IsogeometricApplication_AddCustomUtilities4ToPython()
{

    IsogeometricApplication_AddMeshToPython<2>();
    IsogeometricApplication_AddMeshToPython<3>();

    IsogeometricApplication_AddModelPartToPython<2>();
    IsogeometricApplication_AddModelPartToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos


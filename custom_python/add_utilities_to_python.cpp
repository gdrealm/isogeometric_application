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
//   Date:                $Date: Jan 9, 2013 $
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
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/isogeometric_classical_post_utility.h"
#include "custom_utilities/isogeometric_post_utility.h"
#include "custom_utilities/nurbs_test_utils.h"
#include "custom_utilities/bezier_test_utils.h"
#include "custom_utilities/isogeometric_merge_utility.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_utilities/hierarchical_nurbs/hn_mesh.h"
#include "custom_utilities/domain_manager.h"
#include "custom_utilities/domain_manager_2d.h"

#ifdef ISOGEOMETRIC_USE_HDF5
#include "custom_utilities/hdf5_post_utility.h"
#endif

#ifdef ISOGEOMETRIC_USE_GISMO
#include "custom_utilities/gismo/gismo_mesh.h"
#endif

namespace Kratos
{

namespace Python
{


int BSplineUtils_FindSpan(
    BSplineUtils& dummy,
    const int rN,
    const int rP,
    const double rXi,
    const BSplineUtils::ValuesContainerType& rU
)
{
    return dummy.FindSpan(rN, rP, rXi, rU);
}

void BSplineUtils_BasisFuns(
    BSplineUtils& dummy,
    BSplineUtils::ValuesContainerType& rS,
    const int rI,
    const double rXi,
    const int rP,
    const BSplineUtils::ValuesContainerType& rU
)
{
    dummy.BasisFuns(rS, rI, rXi, rP, rU);
}

void BezierUtils_Bernstein(
    BezierUtils& dummy,
    BezierUtils::ValuesContainerType& rS,
    const int p,
    const double x
)
{
    dummy.bernstein(rS, p, x);
}

void BezierUtils_Bernstein_der(
    BezierUtils& dummy,
    BezierUtils::ValuesContainerType& rS,
    BezierUtils::ValuesContainerType& rD,
    const int p,
    const double x
)
{
    dummy.bernstein(rS, rD, p, x);
}

void BezierUtils_DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(
    BezierUtils& dummy,
    ModelPart::Pointer pModelPart,
    std::string FileName
)
{
    dummy.DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(pModelPart, FileName);
}

template<class T>
void BezierUtils_ComputeCentroid(
    BezierUtils& dummy,
    typename T::Pointer& pElem,
    typename T::GeometryType::PointType& P
)
{
    dummy.ComputeCentroid<T>(pElem, P);
}

void NURBSTestUtils_ProbeGlobalCoordinates2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeGlobalCoordinates(pElement, X, Y);
}

void NURBSTestUtils_ProbeGlobalCoordinates3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeGlobalCoordinates(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeShapeFunctionValues2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeShapeFunctionValues(pElement, X, Y);
}

void NURBSTestUtils_ProbeShapeFunctionValues3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeShapeFunctionValues(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeShapeFunctionDerivatives2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement, X, Y);
}

void NURBSTestUtils_ProbeShapeFunctionDerivatives3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeJacobian2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeJacobian(pElement, X, Y);
}

void NURBSTestUtils_ProbeJacobian3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeJacobian(pElement, X, Y, Z);
}

void IsogeometricApplication_AddCustomUtilitiesToPython()
{
    using namespace boost::python;
    
    enum_<PostElementType>("PostElementType")
    .value("Triangle", _TRIANGLE_)
    .value("Quadrilateral", _QUADRILATERAL_)
    .value("Tetrahedra", _TETRAHEDRA_)
    .value("Hexahedra", _HEXAHEDRA_)
    ;
    
    class_<BSplineUtils, BSplineUtils::Pointer, boost::noncopyable>("BSplineUtils", init<>())
    .def("FindSpan", BSplineUtils_FindSpan)
    .def("BasisFuns", BSplineUtils_BasisFuns)
    .def("test_ComputeBsplinesKnotInsertionCoefficients1DLocal", &BSplineUtils::test_ComputeBsplinesKnotInsertionCoefficients1DLocal)
    ;
    
    class_<BezierUtils, BezierUtils::Pointer, boost::noncopyable>("BezierUtils", init<>())
    .def("Bernstein", BezierUtils_Bernstein)
    .def("BernsteinDerivative", BezierUtils_Bernstein_der)
    .def("DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients", BezierUtils_DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients)
    .def("ComputeCentroid", BezierUtils_ComputeCentroid<Element>)
    .def("ComputeCentroid", BezierUtils_ComputeCentroid<Condition>)
//    .def("compute_extended_knot_vector", &BezierUtils::compute_extended_knot_vector)
//    .def("bezier_extraction_tsplines_1d", &BezierUtils::bezier_extraction_tsplines_1d)
    .def("test_tsplines_1", &BezierUtils::test_tsplines_1)
    .def("test_bezier_extraction_local_1d", &BezierUtils::test_bezier_extraction_local_1d)
    ;
    
    class_<IsogeometricClassicalPostUtility, IsogeometricClassicalPostUtility::Pointer, boost::noncopyable>("IsogeometricClassicalPostUtility", init<ModelPart::Pointer>())
    .def("GenerateModelPart", &IsogeometricClassicalPostUtility::GenerateModelPart)
    .def("GenerateModelPart2", &IsogeometricClassicalPostUtility::GenerateModelPart2)
    .def("GenerateModelPart2AutoCollapse", &IsogeometricClassicalPostUtility::GenerateModelPart2AutoCollapse)
    .def("TransferNodalResults", &IsogeometricClassicalPostUtility::TransferNodalResults<Variable<double> >)
    .def("TransferNodalResults", &IsogeometricClassicalPostUtility::TransferNodalResults<Variable<Vector> >)
    .def("TransferNodalResults", &IsogeometricClassicalPostUtility::TransferNodalResults<Variable<array_1d<double, 3> > >)
    .def("TransferIntegrationPointResults", &IsogeometricClassicalPostUtility::TransferIntegrationPointResults<Variable<double> >)
    .def("TransferIntegrationPointResults", &IsogeometricClassicalPostUtility::TransferIntegrationPointResults<Variable<Vector> >)
    .def("SynchronizeActivation", &IsogeometricClassicalPostUtility::SynchronizeActivation)
    .def("TransferElementalData", &IsogeometricClassicalPostUtility::TransferElementalData<Variable<bool> >)
    .def("TransferConditionalData", &IsogeometricClassicalPostUtility::TransferConditionalData<Variable<bool> >)
    .def("TransferVariablesToNodes", &IsogeometricClassicalPostUtility::TransferVariablesToNodes<Variable<double> >)
    .def("TransferVariablesToNodes", &IsogeometricClassicalPostUtility::TransferVariablesToNodes<Variable<Vector> >)
    .def("GlobalNodalRenumbering", &IsogeometricClassicalPostUtility::GlobalNodalRenumbering)
    ;

    class_<IsogeometricPostUtility, IsogeometricPostUtility::Pointer, boost::noncopyable>("IsogeometricPostUtility", init<>())
    .def("TransferNodalResults", &IsogeometricPostUtility::TransferNodalResults<Variable<double> >)
    .def("TransferNodalResults", &IsogeometricPostUtility::TransferNodalResults<Variable<Vector> >)
    .def("TransferNodalResults", &IsogeometricPostUtility::TransferNodalResults<Variable<array_1d<double, 3> > >)
    .def("TransferIntegrationPointResults", &IsogeometricPostUtility::TransferIntegrationPointResults<Variable<double> >)
    .def("TransferIntegrationPointResults", &IsogeometricPostUtility::TransferIntegrationPointResults<Variable<Vector> >)
    .def("TransferVariablesToNodes", &IsogeometricPostUtility::TransferVariablesToNodes<Variable<double> >)
    .def("TransferVariablesToNodes", &IsogeometricPostUtility::TransferVariablesToNodes<Variable<Vector> >)
    ;

    #ifdef ISOGEOMETRIC_USE_HDF5
    class_<HDF5PostUtility, HDF5PostUtility::Pointer, boost::noncopyable>("HDF5PostUtility", init<const std::string>())
    .def(init<const std::string, const std::string>())
    .def("WriteNodes", &HDF5PostUtility::WriteNodes)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<double>)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<array_1d<double, 3> >)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<Vector>)
    .def("WriteElementalData", &HDF5PostUtility::WriteElementalData<bool>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<double>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<array_1d<double, 3> >)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<Vector>)
    .def("ReadElementalData", &HDF5PostUtility::ReadElementalData<bool>)
    ;
    #endif
    
    class_<NURBSTestUtils, NURBSTestUtils::Pointer, boost::noncopyable>("NURBSTestUtils", init<>())
    .def("Test1", &NURBSTestUtils::Test1)
    .def("Test2", &NURBSTestUtils::Test2)
    .def("ProbeGlobalCoordinates", &NURBSTestUtils_ProbeGlobalCoordinates2)
    .def("ProbeGlobalCoordinates", &NURBSTestUtils_ProbeGlobalCoordinates3)
    .def("ProbeShapeFunctionValues", &NURBSTestUtils_ProbeShapeFunctionValues2)
    .def("ProbeShapeFunctionValues", &NURBSTestUtils_ProbeShapeFunctionValues3)
    .def("ProbeShapeFunctionDerivatives", &NURBSTestUtils_ProbeShapeFunctionDerivatives2)
    .def("ProbeShapeFunctionDerivatives", &NURBSTestUtils_ProbeShapeFunctionDerivatives3)
    .def("ProbeJacobian", &NURBSTestUtils_ProbeJacobian2)
    .def("ProbeJacobian", &NURBSTestUtils_ProbeJacobian3)
    .def("DumpNodalValues", &NURBSTestUtils::DumpNodalValues<double>)
    .def("DumpNodalValues", &NURBSTestUtils::DumpNodalValues<array_1d<double, 3> >)
    ;
    
    class_<IsogeometricMergeUtility, IsogeometricMergeUtility::Pointer, boost::noncopyable>(
        "IsogeometricMergeUtility", init<>())
    .def("Add", &IsogeometricMergeUtility::Add)
    .def("Export", &IsogeometricMergeUtility::Export)
    .def("DumpNodalVariablesList", &IsogeometricMergeUtility::DumpNodalVariablesList)
    ;
    
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

    enum_<HN_ECHO_FLAGS>("HN_ECHO_FLAGS")
    .value("ECHO_REFIMENT", ECHO_REFIMENT)
    ;

    class_<HnMesh, HnMesh::Pointer, boost::noncopyable>
    ("HnMesh", init<std::string>())
    .def("SetEchoLevel", &HnMesh::SetEchoLevel)
    .def("ReadMesh", &HnMesh::ReadMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("Refine", &HnMesh::Refine) // use this for debugging only, use RefineNodes and LinearDependencyRefine instead
    .def("RefineNodes", &HnMesh::RefineNodes)
    .def("LinearDependencyRefine", &HnMesh::LinearDependencyRefine)
    .def("BuildMesh", &HnMesh::BuildMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportCellTopology", &HnMesh::ExportCellTopology)
    .def("ExportCellGeology", &HnMesh::ExportCellGeology)
//    .def("ExportRefinedDomain", &HnMesh::ExportRefinedDomain)
    .def("ExportSupportDomain", &HnMesh::ExportSupportDomain)
    .def("ExportMatlab", &HnMesh::ExportMatlab)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportMDPA", &HnMesh::ExportMDPA)
    .def("ExportMDPA2", &HnMesh::ExportMDPA2)
    .def("ExportPostMDPA", &HnMesh::ExportPostMDPA)
    .def("ExportCellGeologyAsPostMDPA", &HnMesh::ExportCellGeologyAsPostMDPA)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("PrintKnotVectors", &HnMesh::PrintKnotVectors)
    .def("PrintCells", &HnMesh::PrintCells)
    .def("PrintBasisFuncs", &HnMesh::PrintBasisFuncs)
    .def("PrintRefinementHistory", &HnMesh::PrintRefinementHistory)
    .def("CheckNestedSpace", &HnMesh::CheckNestedSpace)
    .def(self_ns::str(self))
    ;

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

    #ifdef ISOGEOMETRIC_USE_GISMO
    class_<GismoMesh, GismoMesh::Pointer, boost::noncopyable>
    ("GismoMesh", init<std::string>())
    .def("SetEchoLevel", &GismoMesh::SetEchoLevel)
    .def("ReadMesh", &GismoMesh::ReadMesh)
    .def(self_ns::str(self))
    ;
    #endif

}

}  // namespace Python.

} // Namespace Kratos


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_MODEL_PART_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_MODEL_PART_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/patch.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

/**
Coupling between KRATOS model_part and multipatch structure
 */
template<int TDim>
class MultiPatchModelPart
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchModelPart);

    /// Type definition
    typedef Element::NodeType NodeType;
    typedef IsogeometricGeometry<NodeType> IsogeometricGeometryType;
    typedef typename Patch<TDim>::ControlPointType ControlPointType;

    /// Default constructor
    MultiPatchModelPart(typename MultiPatch<TDim>::Pointer pMultiPatch)
    : mpMultiPatch(pMultiPatch), mIsModelPartReady(false)
    {
        mpModelPart = ModelPart::Pointer(new ModelPart("MultiPatch"));
    }

    /// Destructor
    virtual ~MultiPatchModelPart() {}

    /// Get the underlying model_part pointer
    ModelPart::Pointer pModelPart() {return mpModelPart;}

    /// Get the underlying model_part pointer
    ModelPart::ConstPointer pModelPart() const {return mpModelPart;}

    /// Get the underlying multipatch pointer
    typename MultiPatch<TDim>::Pointer pMultiPatch() {return mpModelPart;}

    /// Get the underlying multipatch pointer
    typename MultiPatch<TDim>::ConstPointer pMultiPatch() const {return mpModelPart;}

    /// Check if the multipatch model_part ready for transferring/transmitting data
    bool IsReady() const {return mpMultiPatch->IsEnumerated() && mIsModelPartReady;}

    /// Start the process to cook new model_part. This function will first create the new model_part instance and add in the nodes (which are the control points in the multipatch)
    void BeginModelPart()
    {
        mIsModelPartReady = false;

        // check if the multipatch is enumerated, if not enumerate it.
        if (!mpMultiPatch->IsEnumerated())
        {
            mpMultiPatch->Enumerate();
        }

        // create new model_part
        ModelPart::Pointer pNewModelPart = ModelPart::Pointer(new ModelPart(mpModelPart->Name()));

        // create new nodes from control points
        for (std::size_t i = 0; i < mpMultiPatch->EquationSystemSize(); ++i)
        {
            std::tuple<std::size_t, std::size_t> loc = mpMultiPatch->EquationIdLocation(i);

            const std::size_t& patch_id = std::get<0>(loc);
            const std::size_t& local_id = std::get<1>(loc);
//            KRATOS_WATCH(patch_id)
//            KRATOS_WATCH(local_id)

            typedef typename Patch<TDim>::ControlPointType ControlPointType;
            const ControlPointType& point = mpMultiPatch->pGetPatch(patch_id)->pControlPointGridFunction()->pControlGrid()->GetData(local_id);
            // KRATOS_WATCH(point)

            ModelPart::NodeType::Pointer pNewNode = pNewModelPart->CreateNewNode(i+1, point.X(), point.Y(), point.Z());

            // TODO assign corresponding data to new node
        }

        // swap the internal model_part with new model_part
        mpModelPart.swap(pNewModelPart);
    }

    /// create the elements out from the patch and add to the model_part
    void AddElement(typename Patch<TDim>::Pointer pPatch, const std::string& element_name,
            const std::size_t& starting_id, const std::size_t& prop_id)
    {
        if (IsReady()) return;

        // get the grid function for control points
        const GridFunction<TDim, ControlPointType>& rControlPointGridFunction = pPatch->ControlPointGridFunction();

        // construct the cell manager out from the FESpace
        typedef typename FESpace<TDim>::cell_container_t cell_container_t;
        typename cell_container_t::Pointer pCellManager = pPatch->pFESpace()->ConstructCellManager();

        // get the Properties
        Properties::Pointer p_temp_properties = mpModelPart->pGetProperties(prop_id);

        // get the sample element
        if(!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);

        // loop through each cell in the space
        Element::NodesArrayType temp_element_nodes;
        typename IsogeometricGeometryType::Pointer p_temp_geometry;
        std::size_t cnt = starting_id;
        Vector dummy;
        int max_integration_method = 1;
        if (p_temp_properties->Has(NUM_IGA_INTEGRATION_METHOD))
            max_integration_method = (*p_temp_properties)[NUM_IGA_INTEGRATION_METHOD];
        for (typename cell_container_t::iterator it_cell = pCellManager->begin(); it_cell != pCellManager->end(); ++it_cell)
        {
            // KRATOS_WATCH(*(*it_cell))
            // get new nodes
            temp_element_nodes.clear();

            const std::vector<std::size_t>& anchors = (*it_cell)->GetSupportedAnchors();
            Vector weights(anchors.size());
            for (std::size_t i = 0; i < anchors.size(); ++i)
            {
                temp_element_nodes.push_back(( *(FindKey(mpModelPart->Nodes(), anchors[i]+1, "Node").base())));
                weights[i] = rControlPointGridFunction.pControlGrid()->GetData(pPatch->pFESpace()->LocalId(anchors[i])).W();
            }

            // create the geometry
            p_temp_geometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(r_clone_element.GetGeometry().Create(temp_element_nodes));

            p_temp_geometry->AssignGeometryData(dummy,
                                                dummy,
                                                dummy,
                                                weights,
                                                (*it_cell)->GetExtractionOperator(),
                                                static_cast<int>(pPatch->Order(0)),
                                                static_cast<int>(pPatch->Order(1)),
                                                static_cast<int>(pPatch->Order(2)),
                                                max_integration_method);

            // create the element and add to the model_part
            Element::Pointer pNewElement = r_clone_element.Create(cnt++, p_temp_geometry, p_temp_properties);
            mpModelPart->Elements().push_back(pNewElement);
        }

        // sort the element contain and make it consistent
        mpModelPart->Elements().Unique();
    }

    /// create the conditions out from the boundary of the patch and add to the model_part
    void AddCondition(const std::size_t& patch_id, const BoundarySide& side,
            const std::string& condition_name, const std::size_t& starting_id, const std::size_t& prop_id)
    {
        if (IsReady()) return;
    }

    /// Finalize the model_part creation process
    void EndModelPart()
    {
        if (IsReady()) return;
        mIsModelPartReady = true;
    }

    /// Synchronize from multipatch to model_part
    template<class TVariableType>
    void SynchronizeForward()
    {
        if (!IsReady())
            return;
    }

    /// Synchronize from model_part to the multipatch
    template<class TVariableType>
    void SynchronizeBackward()
    {
        if (!IsReady())
            return;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchModelPart";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ">>>ModelPart:" << std::endl;
        rOStream << *mpModelPart << std::endl;
        rOStream << ">>>MultiPatch" << std::endl;
        rOStream << *mpMultiPatch << std::endl;
    }

private:

    bool mIsModelPartReady;

    ModelPart::Pointer mpModelPart;
    typename MultiPatch<TDim>::Pointer mpMultiPatch;

    template<class TContainerType, class TKeyType>
    typename TContainerType::iterator FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName)
    {
        typename TContainerType::iterator i_result;
        if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            std::stringstream buffer;
            buffer << ComponentName << " #" << ThisKey << " is not found.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }

};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchModelPart<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_MODEL_PART_H_INCLUDED


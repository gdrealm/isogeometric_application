//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class MultiPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchUtility);

    /// Type definition

    /// Default constructor
    MultiPatchUtility() {}

    /// Destructor
    virtual ~MultiPatchUtility() {}

    /// Create new patch and wrap it with pointer
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchPointer(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace));
    }

    /// Get the last node id of the model part
    static std::size_t GetLastNodeId(ModelPart& r_model_part)
    {
        std::size_t lastNodeId = 0;
        for(typename ModelPart::NodesContainerType::iterator it = r_model_part.Nodes().begin();
                it != r_model_part.Nodes().end(); ++it)
        {
            if(it->Id() > lastNodeId)
                lastNodeId = it->Id();
        }

        return lastNodeId;
    }


    /// Get the last element id of the model part
    static std::size_t GetLastElementId(ModelPart& r_model_part)
    {
        std::size_t lastElementId = 0;
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
                it != r_model_part.Elements().ptr_end(); ++it)
        {
            if((*it)->Id() > lastElementId)
                lastElementId = (*it)->Id();
        }

        return lastElementId;
    }


    /// Get the last condition id of the model part
    static std::size_t GetLastConditionId(ModelPart& r_model_part)
    {
        std::size_t lastCondId = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            if((*it)->Id() > lastCondId)
                lastCondId = (*it)->Id();
        }

        return lastCondId;
    }


    /// Get the last properties id of the model_part
    static std::size_t GetLastPropertiesId(ModelPart& r_model_part)
    {
        std::size_t lastPropId = 0;
        for(typename ModelPart::PropertiesContainerType::ptr_iterator it = r_model_part.rProperties().ptr_begin();
                it != r_model_part.rProperties().ptr_end(); ++it)
        {
            if((*it)->Id() > lastPropId)
                lastPropId = (*it)->Id();
        }

        return lastPropId;
    }

    /// Find the element in the KRATOS container with specific key
    template<class TContainerType, class TKeyType>
    static typename TContainerType::iterator FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName)
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

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

}; // end class MultiPatchUtility


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED defined


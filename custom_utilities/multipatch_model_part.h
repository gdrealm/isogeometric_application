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
#include "custom_utilities/nurbs/multipatch.h"

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

    /// Default constructor
    MultiPatchModelPart(ModelPart::Pointer mModelPart, typename MultiPatch<TDim>::Pointer pMultiPatch)
    : mpModelPart(pModelPart), mpMultiPatch(pMultiPatch)
    {}

    /// Get the underlying model_part pointer
    ModelPart::Pointer pModelPart() {return mpModelPart;}

    /// Get the underlying model_part pointer
    ModelPart::ConstPointer pModelPart() const {return mpModelPart;}

    /// Get the underlying multipatch pointer
    typename MultiPatch<TDim>::Pointer pMultiPatch() {return mpModelPart;}

    /// Get the underlying multipatch pointer
    typename MultiPatch<TDim>::ConstPointer pMultiPatch() const {return mpModelPart;}

    /// Destructor
    virtual ~MultiPatchModelPart() {}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchModelPart";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    ModelPart::Pointer mpModelPart;
    typename MultiPatch<TDim>::Pointer mpMultiPatch;

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchModelPart& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_MODEL_PART_H_INCLUDED


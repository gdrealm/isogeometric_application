//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/patch.h"

namespace Kratos
{

/**
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<std::size_t TDim>
class MultiPatch
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);

    /// Type definition
    typedef Patch<TDim> PatchType;
    std::vector<typename PatchType::Pointer> PatchContainerType;

    /// Default constructor
    MultiPatch() : mGridSystemSize(0) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPacth(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
    }

    /// Validate the MultiPatch
    virtual bool Validate() const
    {
        for (PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            bool check = (*it)->Validate();
            if (!check)
                return false;
        }

        return true;
    }

    /// Get the grid system size
    const std::size_t& GridSystemSize() const {return mGridSystemSize;}

    /// iterators
    PatchContainerType::iterator begin() {return mpPatches.begin();}
    PatchContainerType::const_iterator begin() const {return mpPatches.begin();}
    PatchContainerType::iterator end() {return mpPatches.end();}
    PatchContainerType::const_iterator end() const {return mpPatches.end();}

    /// Get the number of patches
    const std::size_t& size() const {return mpPatches.size();}

    /// Enumerate all the patches
    void Enumerate(const std::size_t& start)
    {
        // firstly assign all the index to all the patch/grid functions to -1
        for (PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            (*it)->ResetDofs();
        }

        // secondly enumerate each patch
        mGridSystemSize = 0;
        std::set<std::size_t> enumerated_patches;
        for (PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            // enumerate the patch and remember
            mGridSystemSize = (*it)->Enumerate(mGridSystemSize);
            enumerated_patches.insert((*it)->Id());

            // for each patch, enumerate its neighbors
            for (typename PatchType::NeighborPatchContainerType:iterator it2 = (*it)->GetNeighbors().begin();
                    it2 != (*it)->GetNeighbors().end(); ++it2)
            {
                if ( std::find(enumerated_patches.begin(), enumerated_patches.end(), (*it2)->Id()) == enumerated_patches.end() )
                {
                    (*it2)->Enumerate(mGridSystemSize);
                    enumerated_patches.insert((*it2)->Id());
                }
            }
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Number of patches: " << mpPatches.size();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        for (PatchContainerType::iterator it = begin(); it != end(); ++it)
            rOStream << *(*it) << std::endl;
    }

private:

    std::size_t mGridSystemSize;
    PatchContainerType mpPatches;

};

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatch<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED defined


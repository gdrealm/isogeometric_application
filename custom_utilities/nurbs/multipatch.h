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
#include "utilities/indexed_object.h"
#include "containers/pointer_vector_set.h"

namespace Kratos
{

/**
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical NURBS patch, or a T-Splines patch.
 */
template<std::size_t TDim>
class MultiPatch : public boost::enable_shared_from_this<MultiPatch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);

    /// Type definition
    typedef Patch<TDim> PatchType;
    typedef PointerVectorSet<PatchType, IndexedObject> PatchContainerType;

    /// Default constructor
    MultiPatch() : mGridSystemSize(0) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPatch(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
    }

    /// Reset Id for all the patches
    void ResetId()
    {
        std::size_t Id = 0;
        for (typename PatchContainerType::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->SetId(++Id);
        }
    }

    /// Validate the MultiPatch
    virtual bool Validate() const
    {
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            bool check = it->Validate();
            if (!check)
                return false;
        }

        return true;
    }

    /// Get the grid system size
    const std::size_t& GridSystemSize() const {return mGridSystemSize;}

    /// iterators
    typename PatchContainerType::iterator begin() {return mpPatches.begin();}
    typename PatchContainerType::const_iterator begin() const {return mpPatches.begin();}
    typename PatchContainerType::iterator end() {return mpPatches.end();}
    typename PatchContainerType::const_iterator end() const {return mpPatches.end();}

    /// Get the number of patches
    const std::size_t& size() const {return mpPatches.size();}

    /// Enumerate all the patches
    void Enumerate(const std::size_t& start)
    {
        // firstly assign all the index to all the patch/grid functions to -1
        for (typename PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            (*it)->ResetDofs();
        }

        // secondly enumerate each patch
        mGridSystemSize = 0;
        std::set<std::size_t> enumerated_patches;
        for (typename PatchContainerType::iterator it = begin(); it != end(); ++it)
        {
            // enumerate the patch and remember
            mGridSystemSize = (*it)->Enumerate(mGridSystemSize);
            enumerated_patches.insert((*it)->Id());

            // for each patch, enumerate its neighbors
            for (typename PatchType::NeighborPatchContainerType::iterator it2 = (*it)->GetNeighbors().begin();
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

    /// Make the two patches neighbor. This required that two patches are conformed at the interface.
    /// If two patches are conformed, then the grid function on side1 will be transferred to side2.
    static void MakeNeighbor(typename Patch<TDim>::Pointer& pPatch1, const BoundarySide& side1,
            typename Patch<TDim>::Pointer& pPatch2, const BoundarySide& side2)
    {
        typename Patch<TDim-1>::Pointer pBPatch1 = pPatch1->ConstructBoundaryPatch(side1);
        typename Patch<TDim-1>::Pointer pBPatch2 = pPatch2->ConstructBoundaryPatch(side2);

        if( (*pBPatch1) == (*pBPatch2) )
        {
            // synchronize grid function data
            // temporarily disable. I think it can create potential conflict.
            // pBPatch1->SynchronizeGridFunction(*pBPatch2);

            // set the neighbor information
            if (side1 == _LEFT_)
                pPatch1->pSetLeft(pPatch2);
            else if (side1 == _RIGHT_)
                pPatch1->pSetRight(pPatch2);
            else if (side1 == _TOP_)
                pPatch1->pSetTop(pPatch2);
            else if (side1 == _BOTTOM_)
                pPatch1->pSetBottom(pPatch2);
            else if (side1 == _FRONT_)
                pPatch1->pSetFront(pPatch2);
            else if (side1 == _BACK_)
                pPatch1->pSetBack(pPatch2);

            if (side2 == _LEFT_)
                pPatch2->pSetLeft(pPatch1);
            else if (side2 == _RIGHT_)
                pPatch2->pSetRight(pPatch1);
            else if (side2 == _TOP_)
                pPatch2->pSetTop(pPatch1);
            else if (side2 == _BOTTOM_)
                pPatch2->pSetBottom(pPatch1);
            else if (side2 == _FRONT_)
                pPatch2->pSetFront(pPatch1);
            else if (side2 == _BACK_)
                pPatch2->pSetBack(pPatch1);

            // KRATOS_WATCH(*pPatch1)
            // KRATOS_WATCH(*pPatch2)
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The two patch's boundaries are not conformed", "")
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Number of patches: " << mpPatches.size();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
            rOStream << (*it) << std::endl;
    }

protected:

    /// Access the underlying list of patches
    PatchContainerType& Patches() {return mpPatches;}

private:

    std::size_t mGridSystemSize;
    PatchContainerType mpPatches; // container for all the patches

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


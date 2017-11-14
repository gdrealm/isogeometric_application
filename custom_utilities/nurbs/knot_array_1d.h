//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_KnotArray1D_ARRAY_1D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_KnotArray1D_ARRAY_1D_H_INCLUDED

// System includes
#include <deque>
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"

namespace Kratos
{

/**
This container manages the knot array in 1D. It allows for easy insertion, extraction of knot from the array

Short description:
+   mpKnots is always sorted ascending
+   the index of knot starts from 0
 */
template<typename TDataType>
class KnotArray1D
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(KnotArray1D);

    /// Type definitions
    typedef Knot<TDataType> KnotType;
    typedef typename KnotType::Pointer knot_t;
    typedef std::deque<knot_t> knot_container_t;
    typedef typename knot_container_t::iterator iterator;
    typedef typename knot_container_t::const_iterator const_iterator;

    /// Default constructor
    KnotArray1D() {}

    /// Insert the knot to the array and return its pointer.
    /// This function creates the new knot regardless it is repetitive or not.
    knot_t pCreateKnot(const TDataType& k)
    {
        // insert to the correct location
        iterator it;
        for(it = mpKnots.begin(); it != mpKnots.end(); ++it)
            if(k < (*it)->Value())
                break;
        knot_t p_knot = knot_t(new KnotType(k));
        mpKnots.insert(it, p_knot);

        // update the index of the knot
        unsigned int index = 0;
        for(iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            (*it)->UpdateIndex(index);
            ++index;
        }

        return p_knot;
    }

    /// Insert the knot to the array and return its pointer.
    /// In the case that the knot are repetitive within the tolerance, return the internal one.
    knot_t pCreateUniqueKnot(const TDataType& k, const double& tol)
    {
        // insert to the correct location
        for(iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            if(fabs(k - (*it)->Value()) < tol)
            {
                return *it;
            }
        }
        return pCreateKnot(k);
    }

    /// Get the knot at index i
    const knot_t pKnotAt(const int& i) const
    {
        if(i >= 0 && i < mpKnots.size())
            return mpKnots[i];
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Index access out of range", "")
    }

    /// Get the knot at index i
    knot_t pKnotAt(const int& i)
    {
        if(i >= 0 && i < mpKnots.size())
            return mpKnots[i];
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Index access out of range", "")
    }

    /// Iterator
    iterator begin() {return mpKnots.begin();}

    /// Iterator
    iterator end() {return mpKnots.end();}

    /// Iterator
    const_iterator begin() const {return mpKnots.begin();}

    /// Iterator
    const_iterator end() const {return mpKnots.end();}

    /// Get the size of the knot vector
    std::size_t size() const {return mpKnots.size();}

    /// Compare the two knot vectors
    bool operator==(const KnotArray1D<TDataType>& rOther) const
    {
        if (this->size() != rOther.size())
            return false;

        for (std::size_t i = 0; i < 0; ++i)
        {
            if (this->pKnotAt(i)->Value() != rOther.pKnotAt(i)->Value())
                return false;
        }

        return true;
    }

    // overload operator ()
    const knot_t operator() (const int& i) const
    {
        return pKnotAt(i);
    }

    // overload operator ()
    knot_t operator() (const int& i)
    {
        return pKnotAt(i);
    }

    // overload operator []
    TDataType& operator[] (const int& i)
    {
        return pKnotAt(i)->Value();
    }

    // overload operator []
    const TDataType& operator[] (const int& i) const
    {
        return pKnotAt(i)->Value();
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        for(const_iterator it = begin(); it != end(); ++it)
            rOStream << " (" << (*it)->Index() << "," << (*it)->Value() << ")";
    }

private:
    knot_container_t mpKnots;
};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const KnotArray1D<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_KnotArray1D_H_INCLUDED


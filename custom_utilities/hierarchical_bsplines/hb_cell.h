//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 May 2015 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/nurbs/cell.h"

namespace Kratos
{

// forward declaration
class HBBasisFunction;

/**
    Represent a cell in hierarchical B-Splines topology
 */
class HBCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBCell);

    /// Type definitions
    typedef Knot<double>::Pointer knot_t;

    typedef Cell BaseType;

    typedef boost::shared_ptr<HBBasisFunction> bf_t;
    typedef std::set<bf_t> bf_container_t;
    typedef bf_container_t::iterator bf_iterator;

    /// Default constructor
    HBCell(std::size_t Id, unsigned int Level, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp)
    : BaseType(Id, pLeft, pRight, pDown, pUp), mLevel(Level)
    {}

    HBCell(std::size_t Id, unsigned int Level, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp, knot_t pBelow, knot_t pAbove)
    : BaseType(Id, pLeft, pRight, pDown, pUp, pBelow, pAbove), mLevel(Level)
    {}

    /// Destructor
    virtual ~HBCell()
    {}

    /// Get the level of this cell
    unsigned int Level() const {return mLevel;}

    /// Add the basis function to the set. If it does exist in the set, return the internal one
    bf_t AddBf(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            if(*it == p_bf)
                return *it;
        mpBasisFuncs.insert(p_bf);
        return p_bf;
    }

    /// Remove basis function from the set
    void RemoveBf(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            if(*it == p_bf)
            {
                mpBasisFuncs.erase(it);
                break;
            }
    }

    /// Iterators
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_iterator bf_begin() const {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_iterator bf_end() const {return mpBasisFuncs.end();}

    /// Number of bfs contain this cell in the support
    std::size_t size() const {return mpBasisFuncs.size();}

    /// Check if this cell is coverred by another cell
    bool IsCoverred(const HBCell::Pointer p_cell, const int Dim) const
    {
        return BaseType::IsCoverred(BaseType::Pointer(p_cell), Dim);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const;

    virtual void PrintData(std::ostream& rOStream) const;

private:
    unsigned int mLevel;
    bf_container_t mpBasisFuncs; // list of basis functions contain this cell in its support
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HBCell& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED


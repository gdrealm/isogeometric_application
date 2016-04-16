//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HN_CELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HN_CELL_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes 
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "custom_utilities/knot.h"
#include "custom_utilities/cell.h"

namespace Kratos
{

/**
    Represent a HnCell in hierarchical NURBS topology
 */
class HnCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HnCell);

    /// Type definitions
    typedef Cell BaseType;
    typedef HnCell::Pointer cell_t;
    typedef Knot::Pointer knot_t;
    typedef std::map<int, cell_t> cell_container_t;
    typedef cell_container_t::iterator cell_iterator;
    typedef cell_container_t::const_iterator cell_const_iterator;

    /// Default constructor
    HnCell(unsigned int Id, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp)
    : BaseType(Id, pLeft, pRight, pDown, pUp), mIsActive(false)
    {}

    HnCell(unsigned int Id, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp, knot_t pBelow, knot_t pAbove)
    : BaseType(Id, pLeft, pRight, pDown, pUp, pBelow, pAbove), mIsActive(false)
    {}

    /// Destructor
    virtual ~HnCell() {}

    /// Check if this cell is active
    bool IsActive() const {return mIsActive;}

    /// Activation and deactivation
    void Activate() {mIsActive = true;}
    void Deactivate() {mIsActive = false;}

    /// Add a child cell to this cell; the child is defined as the cells contained in this cell
    void AddChild(cell_t p_cell) {mpCells[p_cell->Id()] = p_cell;}

    /// Get number of children of this cell
    std::size_t NumberOfChildren() const {return mpCells.size();}

    /// Iterators to the children cells
    cell_iterator cell_begin() {return mpCells.begin();}
    cell_iterator cell_end() {return mpCells.end();}
    cell_const_iterator cell_const_begin() const {return mpCells.begin();}
    cell_const_iterator cell_const_end() const {return mpCells.end();}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }

private:
    bool mIsActive;
    cell_container_t mpCells; // list of children cells of this cell in the next level
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HnCell& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " supporting basis functions: ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_CELL_H_INCLUDED


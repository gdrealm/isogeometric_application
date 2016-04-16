//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HN_LEVEL_2D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HN_LEVEL_2D_H_INCLUDED

// System includes
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

// External includes 
#include <omp.h>
#include "boost/progress.hpp"
#include "boost/algorithm/string.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "hn_basis_function.h"
#include "hn_cell.h"
#include "custom_utilities/knot.h"

namespace Kratos
{

/**
    Class represents a Hierarchical NURBS mesh in 2D

*/
class HnLevel
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HnLevel);

    /// Type definition
    typedef Knot::Pointer knot_t;
    typedef HnCell::Pointer cell_t;
    typedef HnBasisFunction::Pointer bf_t;
    typedef std::map<int, cell_t> cell_container_t;
    typedef cell_container_t::iterator cell_iterator;
    typedef cell_container_t::const_iterator cell_const_iterator;
    typedef std::vector<bf_t> bf_container_t;
    typedef bf_container_t::iterator bf_iterator;
    typedef bf_container_t::const_iterator bf_const_iterator;
    typedef std::deque<knot_t> knot_container_t;

    /// Default constructor
    HnLevel(int Id) : mId(Id) {}

    /// Destructor
    ~HnLevel() {}

    /// Get the Id of this level
    unsigned int Id() const {return mId;}

    /// Get the number of basis functions at this level
    std::size_t NumberOfBasisFunctions() const {return mpBasisFuncs.size();}

    /// Get the number of cells at this level
    std::size_t NumberOfCells() const {return mpCells.size();}

    /// Iterators to the basis function of this level
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_const_iterator bf_const_begin() const {return mpBasisFuncs.begin();}
    bf_const_iterator bf_const_end() const {return mpBasisFuncs.end();}

    /// Iterators to the cell of this level
    cell_iterator cell_begin() {return mpCells.begin();}
    cell_iterator cell_end() {return mpCells.end();}
    cell_const_iterator cell_const_begin() const {return mpCells.begin();}
    cell_const_iterator cell_const_end() const {return mpCells.end();}

    /// Print information of this hierarchical mesh
    void PrintInfo(std::ostream& rOStream) const
    {
        int num_active = 0;
        for(bf_const_iterator it = bf_const_begin(); it != bf_const_end(); ++it)
            if((*it)->IsActive())
                ++num_active;

        rOStream << "Level " << Id() << ": "
                 << mpCells.size() << " cell(s), "
                 << mpBasisFuncs.size() << " basis function(s), "
                 << num_active << " active basis function(s)";
    }

    /// Print data of this hierarchical mesh
    void PrintData(std::ostream& rOStream) const
    {
        // Print the basis functions
        rOStream << "Level " << Id() << " Basis Function list:" << std::endl;
        for(bf_const_iterator it = bf_const_begin(); it != bf_const_end(); ++it)
            rOStream << *(*it) << std::endl;
    }

    /// Initialize the internal data for this level
    /// Remarks: knot vectors must be sorted before calling this function
    void Initialize(unsigned int FirstCellId,
                    unsigned int FirstBasisFuncId,
                    int p1, int p2,
                    knot_container_t pKnots1,
                    knot_container_t pKnots2)
    {
        // firstly compute the number of univariate basis function in each direction
        int num1 = pKnots1.size() - p1 - 1;
        int num2 = pKnots2.size() - p2 - 1;
        
//        std::cout << "pKnots1:";
//        for(int i = 0; i < pKnots1.size(); ++i)
//            std::cout << " " << pKnots1[i]->Value();
//        std::cout << std::endl;
//        std::cout << "pKnots2:";
//        for(int i = 0; i < pKnots2.size(); ++i)
//            std::cout << " " << pKnots2[i]->Value();
//        std::cout << std::endl;
        
        // secondly create list of cells. The cell covers the knot spans
        unsigned int CellId = FirstCellId;
        knot_t pLeft, pRight, pUp, pDown;
        cell_t pCell;
        int cnt = 0;
        double area, tol = 1.0e-6;
        for(std::size_t j = 0; j < pKnots2.size()-1; ++j)
            for(std::size_t i = 0; i < pKnots1.size()-1; ++i)
            {
                pLeft  = pKnots1[i];
                pRight = pKnots1[i + 1];
                pDown  = pKnots2[j];
                pUp    = pKnots2[j + 1];
                // check if the cell is not zero area
                area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                if(fabs(area) > tol)
                {
                    pCell = cell_t(new cell_t::element_type(++CellId, pLeft, pRight, pDown, pUp));
                    mpCells[cnt] = pCell;
                }
                ++cnt;
            }

        // thirdly create the basis function and add supporting cells. Note that the basis function supporting domain is (p+1)x(q+1) cells
        unsigned int BasisFuncId = FirstBasisFuncId;
        HnBasisFunction::Pointer pBasisFunc;
        for(std::size_t j = 0; j < num2; ++j)
            for(std::size_t i = 0; i < num1; ++i)
            {
                // add cells to the basis function
//                int ifunc = j * num1 + i;
                pBasisFunc = HnBasisFunction::Pointer(new HnBasisFunction(++BasisFuncId));
                for(std::size_t l = 0; l < p2 + 1; ++l)
                    for(std::size_t k = 0; k < p1 + 1; ++k)
                    {
                        int icell = (i+k) + (j+l) * (pKnots1.size()-1);
                        if(mpCells.find(icell) != mpCells.end())
                            pBasisFunc->AddCell(mpCells[icell]);
                    }

                // add knots to the basis function
                std::vector<double> knots1;
                std::vector<double> knots2;
                for(std::size_t k = 0; k < p1 + 2; ++k)
                    knots1.push_back(pKnots1[i + k]->Value());
                for(std::size_t k = 0; k < p2 + 2; ++k)
                    knots2.push_back(pKnots2[j + k]->Value());
                pBasisFunc->SetLocalKnotVectors(knots1, knots2);

                // add basis functions to the list
                mpBasisFuncs.push_back(pBasisFunc);
                
//                std::cout << "bf " << pBasisFunc->Id() << " cells:";
//                for(HnBasisFunction::cell_iterator it = pBasisFunc->cell_begin(); it != pBasisFunc->cell_end(); ++it)
//                    std::cout << " " << it->second->Id();
//                std::cout << std::endl;
            }
    }

    /// Initialize the internal data for this level
    /// Remarks: knot vectors must be sorted before calling this function
    void Initialize(int FirstCellId,
                    int FirstBasisFuncId,
                    int p1, int p2, int p3,
                    knot_container_t pKnots1,
                    knot_container_t pKnots2,
                    knot_container_t pKnots3)
    {
        // TODO
    }

    /// Set the coordinates for the basis functions
    void SetCoordinates(int Num1, int Num2,
                        std::vector<double>& Xcoords,
                        std::vector<double>& Ycoords,
                        std::vector<double>& Weights)
    {
        int ifunc;
        for(std::size_t j = 0; j < Num2; ++j)
            for(std::size_t i = 0; i < Num1; ++i)
            {
                ifunc = j * Num1 + i;
                mpBasisFuncs[ifunc]->SetCoordinates(Xcoords[ifunc], Ycoords[ifunc], 0.0, Weights[ifunc]);
            }
    }

    /// Set the coordinates for the basis functions
    void SetCoordinates(int Num1, int Num2, int Num3,
                        std::vector<double>& Xcoords,
                        std::vector<double>& Ycoords,
                        std::vector<double>& Zcoords,
                        std::vector<double>& Weights)
    {
        int ifunc;
        for(std::size_t k = 0; k < Num3; ++k)
            for(std::size_t j = 0; j < Num2; ++j)
                for(std::size_t i = 0; i < Num1; ++i)
                {
                    ifunc = (k * Num2 + j) * Num1 + i;
                    mpBasisFuncs[ifunc]->SetCoordinates(Xcoords[ifunc], Ycoords[ifunc], Zcoords[ifunc], Weights[ifunc]);
                }
    }

private:
    unsigned int mId;
    cell_container_t mpCells;
    bf_container_t mpBasisFuncs;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HnLevel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_LEVEL_2D_H_INCLUDED


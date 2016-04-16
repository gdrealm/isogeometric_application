//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HN_BASIS_FUNCTION_2D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HN_BASIS_FUNCTION_2D_H_INCLUDED

// System includes
#include <ctime>
#include <cmath>
#include <climits>
#include <string>
#include <vector>
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
#include "custom_utilities/bezier_utils.h"
#include "hn_cell.h"

namespace Kratos
{

/**
    Class represents a basis function in Hierarchical NURBS mesh
*/
class HnBasisFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HnBasisFunction);

    /// Type definition
    typedef HnCell::Pointer cell_t;
    typedef std::vector<HnBasisFunction::Pointer> bf_container_t;
    typedef bf_container_t::iterator bf_iterator;
    typedef bf_container_t::const_iterator bf_const_iterator;
    typedef std::vector<cell_t> cell_container_t;
    typedef cell_container_t::iterator cell_iterator;
    typedef cell_container_t::const_iterator cell_const_iterator;

    /// Default constructor
    HnBasisFunction(unsigned int Id) : mId(Id), mIsActive(false), mX(0.0), mY(0.0), mZ(0.0), mW(0.0) {}

    /// Destructor
    ~HnBasisFunction() {}

    /**************************************************************************
                            MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Get the number of children of this basis function
    std::size_t NumberOfChildren() const {return mpChilds.size();}

    /// Add a child which support this basis function
    void AddChild(HnBasisFunction::Pointer pBasisFunc, double RefinedCoefficient)
    {
        mpChilds.push_back(pBasisFunc);
        mRefinedCoefficients[pBasisFunc->Id()] = RefinedCoefficient;
    }

    /// Add a cell support this basis function to the list
    void AddCell(cell_t p_cell) {mpCells.push_back(p_cell);}

    /// Check if this basis function is active or inactive
    bool IsActive() const {return mIsActive;}

    /// Activate this basis function
    void Activate() {mIsActive = true;}

    /// Deactivate this basis function
    void Deactivate() {mIsActive = false;}

    /// Set the local knot vectors to this basis function
    void SetLocalKnotVectors(const std::vector<double>& rKnots1, const std::vector<double>& rKnots2)
    {
        mLocalKnots1.resize(rKnots1.size());
        mLocalKnots2.resize(rKnots2.size());
        std::copy(rKnots1.begin(), rKnots1.end(), mLocalKnots1.begin());
        std::copy(rKnots2.begin(), rKnots2.end(), mLocalKnots2.begin());
    }

    /**************************************************************************
                            ACCESS SUBROUTINES
    **************************************************************************/

    /// Iterators to the child of this basis function
    bf_iterator bf_begin() {return mpChilds.begin();}
    bf_iterator bf_end() {return mpChilds.end();}
    bf_const_iterator bf_const_begin() const {return mpChilds.begin();}
    bf_const_iterator bf_const_end() const {return mpChilds.end();}

    /// Iterators to the supporting cell of this basis function
    cell_iterator cell_begin() {return mpCells.begin();}
    cell_iterator cell_end() {return mpCells.end();}
    cell_const_iterator cell_const_begin() const {return mpCells.begin();}
    cell_const_iterator cell_const_end() const {return mpCells.end();}

    /// Get the Id of this basis function. Each basis function should have unique Id
    unsigned int Id() const {return mId;}

    /// Get and Set the coordinates
    double X0() const {return mX/mW;}
    double Y0() const {return mY/mW;}
    double Z0() const {return mZ/mW;}
    double W() const {return mW;}
    void SetCoordinates(double X, double Y, double Z, double W) {mX = W*X; mY = W*Y; mZ = W*Z; mW = W;}
    void AddCoordinates(double X, double Y, double Z, double W) {mX += W*X; mY += W*Y; mZ += W*Z; mW += W;}

    /// Get the local knot vectors
    template<class ValuesContainerType>
    void GetLocalKnots(int dim, ValuesContainerType& rKnots) const
    {
        if(dim == 1)
        {
            if(rKnots.size() != mLocalKnots1.size())
                rKnots.resize(mLocalKnots1.size());
            std::copy(mLocalKnots1.begin(), mLocalKnots1.end(), rKnots.begin());
        }
        else if(dim == 2)
        {
            if(rKnots.size() != mLocalKnots2.size())
                rKnots.resize(mLocalKnots2.size());
            std::copy(mLocalKnots2.begin(), mLocalKnots2.end(), rKnots.begin());
        }
        else if(dim == 3)
        {
            if(rKnots.size() != mLocalKnots3.size())
                rKnots.resize(mLocalKnots3.size());
            std::copy(mLocalKnots3.begin(), mLocalKnots3.end(), rKnots.begin());
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension parameters", "")
    }

    /// Get the bounding box (=support domain) of this basis function
    void GetBoundingBox(double& Xmin, double& Xmax, double& Ymin, double& Ymax, double& Zmin, double& Zmax)
    {
        Xmin = static_cast<double>(INT_MAX);
        Xmax = -Xmin;
        Ymin = Xmin;
        Ymax = -Ymin;
        Zmin = Xmin;
        Zmax = -Zmin;
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if((*it)->LeftValue() < Xmin)
                Xmin = (*it)->LeftValue();
            if((*it)->RightValue() > Xmax)
                Xmax = (*it)->RightValue();
            if((*it)->DownValue() < Ymin)
                Ymin = (*it)->DownValue();
            if((*it)->UpValue() > Ymax)
                Ymax = (*it)->UpValue();
            if((*it)->BelowValue() < Zmin)
                Zmin = (*it)->BelowValue();
            if((*it)->AboveValue() > Zmax)
                Zmax = (*it)->AboveValue();
        }
    }

    /// Get the refined coefficient of a child
    double GetRefinedCoefficient(int ChildId) const
    {
        std::map<int, double>::const_iterator it = mRefinedCoefficients.find(ChildId);
        if(it != mRefinedCoefficients.end())
            return it->second;
        else
        {
            std::stringstream ss;
            ss << "The basis function " << ChildId << " is not the child of basis function " << Id();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /**************************************************************************
                            COMPUTATION SUBROUTINES
    **************************************************************************/

    /// Compute the Bezier extraction operator of this basis function on the cell
    void ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2);
    void ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2, int p3);

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Print information of this basis function
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Bf(id:" << Id()
                 << ",(x:" << X0()
                 << ",y:"  << Y0()
                 << ",z:"  << Z0()
                 << ",w:"  << W()
                 << "))";
    }

    /// Print data of this basis function
    void PrintData(std::ostream& rOStream) const
    {
        // Print the cells
        rOStream << "Supporting cells:" << std::endl;
        for(cell_const_iterator it = cell_const_begin(); it != cell_const_end(); ++it)
            rOStream << *(*it) << std::endl;
        rOStream << "List of children:";
        for(bf_const_iterator it = bf_const_begin(); it != bf_const_end(); ++it)
        {
            std::map<int, double>::const_iterator it_coeff = mRefinedCoefficients.find((*it)->Id());
            rOStream << " (" << (*it)->Id() << "," << it_coeff->second << ")";
        }
        if(bf_const_end() == bf_const_begin())
            rOStream << " none";
        rOStream << std::endl;
    }

private:
    unsigned int mId;
    bool mIsActive;
    double mX, mY, mZ, mW; // homogeneous coordinates
    bf_container_t mpChilds; // list of refined basis functions that constitute this basis function
    std::map<int, double> mRefinedCoefficients; // coefficient of refined basis functions
    cell_container_t mpCells; // list of cells support this basis function at the level of this basis function
    std::vector<double> mLocalKnots1;
    std::vector<double> mLocalKnots2;
    std::vector<double> mLocalKnots3;
    
    /// Find the span of knot in the local knot vector
    /// Remarks: it will give the based-1 index
    ///          it only works if U is non-repeated
    template<class ValuesContainerType>
    int FindSpanLocal(double Xi, const ValuesContainerType& U)
    {
        if(!U.empty())
        {
            if(Xi < U[0])
                return 0;

            if(Xi > U[U.size()-1])
                return U.size();
            
            for(std::size_t i = 0; i < U.size()-1; ++i)
                if(Xi >= U[i] and Xi < U[i + 1])
                    return i + 1;
        }
        
        return 0;
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HnBasisFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_BASIS_FUNCTION_2D_H_INCLUDED


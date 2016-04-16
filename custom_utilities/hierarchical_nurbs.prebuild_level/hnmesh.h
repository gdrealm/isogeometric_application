//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HN_MESH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HN_MESH_H_INCLUDED

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
#include "boost/foreach.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/knot.h"
#include "hn_basis_function.h"
#include "hn_level.h"

namespace Kratos
{

/**
    Class represents a Hierarchical NURBS mesh in 2D

*/
class HnMesh
{
private:
    template<class data_type>
    class unordered_pair
    {
        public:
            unordered_pair(data_type i, data_type j)
            {
                if(i < j)
                {
                    mi = i;
                    mj = j;
                }
                else
                {
                    mi = j;
                    mj = i;
                }
            }
            
            ~unordered_pair() {}
            
            bool operator<(const unordered_pair& rOther) const
            {
                if(mi == rOther.mi)
                    return mj < rOther.mj;
                else
                    return mi < rOther.mi;
            }
            
            data_type first() const {return mi;}
            data_type second() const {return mj;}

        private:
            data_type mi, mj;
    };

    template<class data_type>
    class unordered_tuple
    {
        public:
            unordered_tuple(data_type i, data_type j, data_type k, data_type l)
            {
                // TODO
            }
            
            ~unordered_tuple() {}
            
            bool operator<(const unordered_tuple& rOther) const
            {
                // TODO
            }
            
        private:
            data_type mi, mj, mk, ml;
    };

public:
    /// Const definition
    static const int NO_READ          = 0;
    static const int READ_PATCH       = 1;
    static const int READ_ORDER       = 2;
    static const int READ_NUMBER      = 3;
    static const int READ_KNOTS       = 4;
    static const int READ_COORDINATES = 5;
    static const int READ_WEIGHTS     = 6;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HnMesh);

    /// Type definition
    typedef Knot<double>::Pointer knot_t;
    typedef std::deque<knot_t> knot_container_t;
    typedef HnLevel::Pointer level_t;
    typedef HnBasisFunction::Pointer bf_t;
    typedef HnCell::Pointer cell_t;
    typedef std::vector<level_t> level_container_t;
    typedef level_container_t::iterator level_iterator;
    typedef level_container_t::reverse_iterator level_reverse_iterator;
    typedef level_container_t::const_iterator level_const_iterator;
    typedef level_container_t::const_reverse_iterator level_const_reverse_iterator;

    /// Default constructor
    HnMesh(std::string Name) : mName(Name), mDebugLevel(1) {}

    /// Destructor
    ~HnMesh() {}

    /// Get the name of this hierarchical mesh
    std::string Name() const {return mName;}

    /// Set the debug level for this mesh
    void SetDebugLevel(int Level) {mDebugLevel = Level;}

    /// Print information of this hierarchical mesh
    void PrintInfo(std::ostream& rOStream) const;

    /// Print data of this hierarchical mesh
    void PrintData(std::ostream& rOStream) const;

    /// Add a new level
    void AddLevel(level_t p_level) {mpLevels.push_back(p_level);}

    /// Get the number of levels in this mesh
    std::size_t NumberOfLevels() const {return mpLevels.size();}

    /// Iterators to the levels of the hierarchical mesh
    level_iterator level_begin() {return mpLevels.begin();}
    level_iterator level_end() {return mpLevels.end();}
    level_reverse_iterator level_reverse_begin() {return mpLevels.rbegin();}
    level_reverse_iterator level_reverse_end() {return mpLevels.rend();}
    level_const_iterator level_const_begin() const {return mpLevels.begin();}
    level_const_iterator level_const_end() const {return mpLevels.end();}
    level_const_reverse_iterator level_const_reverse_begin() const {return mpLevels.rbegin();}
    level_const_reverse_iterator level_const_reverse_end() const {return mpLevels.rend();}

    /// Read and construct the first level from file
    void ReadMesh(std::string fn);

    /// Create the next level and construct the basis function linkages
    void BuildLevel();

    /// Refine on a cuboid domain
    void RefineWindow(double Xi_min, double Xi_max, double Eta_min, double Eta_max, double Zeta_min, double Zeta_max);

    /// Refine the basis function on the nodes given by Id
    void RefineNodes(boost::python::list& pyList);

    /// Refine compatibly (ensure linear independence) the basis functions in all level
    void Refine(std::map<unsigned int, std::vector<bf_t> >& ToBeRefinedBfs);

    /// Export the cell topology to Matlab for rendering
    void ExportCellTopology(std::string fn, bool cell_numbering) const;

    /// Export the basis function topology to Matlab for rendering
    void ExportBasisFunctionTopology(std::string fn, bool bf_numbering, boost::python::list& pyList) const;

    /// Export the cell geological map to Matlab for rendering
    void GenerateCellGeology(std::vector<unsigned int>& point_list, // list of node Id (of the active basis function) in the hierarchical mesh
                             std::map<unsigned int, double>& X_list, // global x-coordinate of the node in point_list
                             std::map<unsigned int, double>& Y_list, // global y-coordinate of the node in point_list
                             std::map<unsigned int, double>& Z_list, // global z-coordinate of the node in point_list
                             std::map<unsigned int, double>& xi_list, // local xi-coordinate of the node in point_list
                             std::map<unsigned int, double>& eta_list, // local eta-coordinate of the node in point_list
                             std::map<unsigned int, double>& zeta_list, // local eta-coordinate of the node in point_list
                             std::map<unsigned int, unsigned int>& cell_list, // cell index of the node in point_list
                             std::set<unordered_pair<unsigned int> >& Lines, // list of lines for post-processing
                             std::set<unordered_tuple<unsigned int> >& Faces, // list of faces for post-processing
                             std::map<unsigned int, std::vector<std::vector<unsigned int> > >& Connectivities // element connectivities
                             ) const;
    void ExportCellGeology(std::string fn) const;
    void ExportCellGeologyAsPostMDPA(std::string fn) const;

    /// Export the information (basis functions, Bezier elements) to Matlab for plotting and debugging
    void ExportMatlab(std::string fn) const;

    /// Export this hierarchical mesh to MDPA
    void ExportMDPA(std::string fn) const;

    /// Export the post model part
    void ExportPostMDPA(std::string fn, int Division1, int Division2, int Division3) const;

    /****************Debugging Utilities*****************/
    /// Find and print all the active basis function in a domain; use for debugging
    void FindBfs(unsigned int level, double xi_min, double xi_max, double eta_min, double eta_max, double zeta_min, double zeta_max) const;

    /// Print all the child of a bf
    void PrintBfChildren(unsigned int Id) const;

    /// Print local knot vector of a bf
    void PrintLocalKnotVector(unsigned int Id) const;

private:
    std::string mName;
    unsigned int mDim;
    unsigned int mOrder1;
    unsigned int mOrder2;
    unsigned int mOrder3;
    int mDebugLevel;
    knot_container_t mpKnots1; // container for knot vector in u-direction
    knot_container_t mpKnots2; // container for knot vector in v-direction
    knot_container_t mpKnots3; // container for knot vector in v-direction
    level_container_t mpLevels; // container for all level of this mesh

    /// ask cell to check if it supports this basis function; otherwise call its children recursively
    void RecursiveAddSupportedBasisFunction(cell_t p_cell, bf_t p_bf);

    /// check if a basis function is linear dependent on a set of sub-basis function in the hierarchical mesh
    /// return false if this bf is not able to be represented by a linear combination of a number of sub-bf in the next levels
    /// return true if it can be represented
    bool IsLinearDependent(bf_t p_bf) const;

    /// check if a cell is coverred by a set of its sub-cells
    bool IsCoverred(cell_t p_cell) const;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HnMesh& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

// #undef DEBUG_LEVEL1
// #undef DEBUG_LEVEL2
// #undef ENABLE_PROFILING

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_MESH_2D_H_INCLUDED


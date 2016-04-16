#include <fstream>
#include <algorithm>
#include "utilities/openmp_utils.h"
#include "hnmesh.h"
// #include "hn_cuboid.h" //unused
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/triangulation_utils.h"
#include "utilities/auto_collapse_spatial_binning.h"

// Warning: currently, enable this will render the isogeometric_post_utility useless
// #define MDPA_NODE_RENUMBERING
// #define MDPA_CELL_RENUMBERING
#define ENABLE_PROFILING
#define DEBUG_REFINE

namespace Kratos
{

    void HnMesh::PrintInfo(std::ostream& rOStream) const
    {
        // Print the knot vector
        rOStream << "knot vector 1:";
        for(knot_container_t::const_iterator it = mpKnots1.begin(); it != mpKnots1.end(); ++it)
            rOStream << " " << (*it)->Value();
        rOStream << std::endl;

        rOStream << "knot vector 2:";
        for(knot_container_t::const_iterator it = mpKnots2.begin(); it != mpKnots2.end(); ++it)
            rOStream << " " << (*it)->Value();
        rOStream << std::endl;

        rOStream << "knot vector 3:";
        for(knot_container_t::const_iterator it = mpKnots3.begin(); it != mpKnots3.end(); ++it)
            rOStream << " " << (*it)->Value();
        rOStream << std::endl;

        // Print the level
        rOStream << "Level summary:" << std::endl;
        for(level_const_iterator it = level_const_begin(); it != level_const_end(); ++it)
        {
            (*it)->PrintInfo(rOStream);
            rOStream << std::endl;
        }
    }

    void HnMesh::PrintData(std::ostream& rOStream) const
    {
        if(mDebugLevel > 2)
        {
            // Print the level
            rOStream << "Level information:" << std::endl;
            for(level_const_iterator it = level_const_begin(); it != level_const_end(); ++it)
            {
                (*it)->PrintData(rOStream);
                rOStream << "------------------------------" << std::endl;;
            }
        }
        
        if(mDebugLevel > 1)
        {
            // Print the active basis functions on each level
            rOStream << "Active basis functions on each level:" << std::endl;
            for(level_const_iterator it = level_const_begin(); it != level_const_end(); ++it)
            {
                rOStream << "Level " << (*it)->Id() << ":" << std::endl;
                int num_active = 0;
                for(HnLevel::bf_const_iterator it2 = (*it)->bf_const_begin(); it2 != (*it)->bf_const_end(); ++it2)
                {
                    if((*it2)->IsActive())
                    {
                        rOStream << *(*it2);
                        rOStream << "--------------" << std::endl;
                        ++num_active;
                    }
                }
                rOStream << "Level " << (*it)->Id() << " has " << num_active << " active basis function(s)" << std::endl;
                rOStream << "------------------------------" << std::endl;
            }
        }
    
    }

    void HnMesh::ReadMesh(std::string fn)
    {
        std::ifstream infile(fn.c_str());
        if(!infile)
            KRATOS_THROW_ERROR(std::logic_error, "Error open file", fn)

        std::string line;
        std::vector<std::string> words;
        int read_mode = READ_PATCH;
        int npatches, dim_index = 0;
        std::vector<int> orders;
        std::vector<int> numbers;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        std::vector<double> z_coords;
        std::vector<double> weights;
        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);
            
            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == READ_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                    }

                    // read info
                    mDim = atoi(words[0].c_str());
                    npatches = atoi(words[1].c_str());
                    if(npatches > 1)
                        KRATOS_THROW_ERROR(std::logic_error, "Currently number of patches > 1 is not supported", "")
                    read_mode = READ_ORDER;
                    continue;
                }
                
                if(read_mode == READ_ORDER)
                {
                    // bound check
                    if(words.size() != mDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Order section must contained number of information equal to dimension, current number of information =", words.size())

                    // read info
                    for(std::size_t i = 0; i < mDim; ++i)
                        orders.push_back(atoi(words[i].c_str()));
                    read_mode = READ_NUMBER;
                    continue;
                }
                
                if(read_mode == READ_NUMBER)
                {
                    // bound check
                    if(words.size() != mDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Number section must contained number of information equal to dimension, current number of information =", words.size())

                    for(std::size_t i = 0; i < mDim; ++i)
                        numbers.push_back(atoi(words[i].c_str()));
                    read_mode = READ_KNOTS;
                    continue;
                }

                if(read_mode == READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers[dim_index] + orders[dim_index] + 1;
                    if(words.size() != knot_len)
                        KRATOS_THROW_ERROR(std::logic_error, "The Knots section must contained number of information equal to n+p+1, current number of information =", words.size())

                    for(std::size_t i = 0; i < knot_len; ++i)
                    {
                        knot_t p_knot = knot_t(new Knot(atof(words[i].c_str())));
                        if(dim_index == 0)
                            mpKnots1.push_back(p_knot);
                        else if(dim_index == 1)
                            mpKnots2.push_back(p_knot);
                        else
                            KRATOS_THROW_ERROR(std::logic_error, "Wrong knot dimension index. Something must be wrong", "")
                    }
                    
                    ++dim_index;
                    if(dim_index == mDim)
                    {
                        dim_index = 0;
                        read_mode = READ_COORDINATES;
                    }
                    continue;
                }
                
                if(read_mode == READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < mDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Coordinates section must contained number of information equal to prod(ni), current number of information =", words.size())

                    if(dim_index == 0)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            x_coords.push_back(atof(words[i].c_str()));
                    else if(dim_index == 1)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            y_coords.push_back(atof(words[i].c_str()));
                    else if(dim_index == 2)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            z_coords.push_back(atof(words[i].c_str()));
                    
                    ++dim_index;
                    if(dim_index == mDim)
                    {
                        dim_index = 0;
                        read_mode = READ_WEIGHTS;
                    }
                    continue;
                }
                
                if(read_mode == READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < mDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Weights section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        weights.push_back(atof(words[i].c_str()));

                    read_mode = NO_READ;
                    continue;
                }
            }
        }

        infile.close();

        // update the order
        if(orders.size() > 0)
            mOrder1 = orders[0];
        if(orders.size() > 1)
            mOrder2 = orders[1];
        if(orders.size() > 2)
            mOrder3 = orders[2];

        // update the index of the knots
        for(std::size_t i = 0; i < mpKnots1.size() ; ++i)
            mpKnots1[i]->UpdateIndex(i);
        for(std::size_t i = 0; i < mpKnots2.size() ; ++i)
            mpKnots2[i]->UpdateIndex(i);
        for(std::size_t i = 0; i < mpKnots3.size() ; ++i)
            mpKnots3[i]->UpdateIndex(i);

        // rescale the coordinates
        for(std::size_t i = 0; i < weights.size(); ++i)
        {
            x_coords[i] /= weights[i];
            y_coords[i] /= weights[i];
            if(mDim > 2)
                z_coords[i] /= weights[i];
        }

        // create and initialize the first level
        level_t pLevel = level_t(new HnLevel(mpLevels.size() + 1));
        if(mDim == 2)
        {
            pLevel->Initialize(0, 0, mOrder1, mOrder2, mpKnots1, mpKnots2);
            pLevel->SetCoordinates(numbers[0], numbers[1], x_coords, y_coords, weights);
        }
        else if(mDim == 3)
        {
            pLevel->Initialize(0, 0, mOrder1, mOrder2, mOrder3, mpKnots1, mpKnots2, mpKnots3);
            pLevel->SetCoordinates(numbers[0], numbers[1], numbers[2], x_coords, y_coords, z_coords, weights);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
        AddLevel(pLevel);

        // activate all basis functions in the first level
        for(HnLevel::bf_iterator it = (*level_begin())->bf_begin(); it != (*level_begin())->bf_end(); ++it)
            (*it)->Activate();

        // activate all cells in the first level
        for(HnLevel::cell_iterator it = (*level_begin())->cell_begin(); it != (*level_begin())->cell_end(); ++it)
            it->second->Activate();

        // for all the cells, check all the supported basis functions to this cell and add to list
        for(HnLevel::bf_iterator it_bf = (*level_begin())->bf_begin(); it_bf != (*level_begin())->bf_end(); ++it_bf)
        {
            if((*it_bf)->IsActive())
            {
                for(HnBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                {
                    // add the supported basis functions and compute the Bezier extraction operator over cells
                    RecursiveAddSupportedBasisFunction(*it_cell, *it_bf);
                }
            }
        }
                
        std::cout << "ReadMesh and build level " << pLevel->Id() << " completed" << std::endl;
    }
    
    void HnMesh::BuildLevel()
    {
        #ifdef ENABLE_PROFILING
        double time_start = OpenMPUtils::GetCurrentTime();
        #endif
    
        // extract all the distinguished knots in each direction
        std::set<double> knot_set1; // set of distinguished knots
        std::vector<double> knots1; // full knot vector
        std::set<double> knot_set2;
        std::vector<double> knots2;
        std::set<double> knot_set3;
        std::vector<double> knots3;
        for(std::size_t i = 0; i < mpKnots1.size(); ++i)
        {
            knot_set1.insert(mpKnots1[i]->Value());
            knots1.push_back(mpKnots1[i]->Value());
        }
        for(std::size_t i = 0; i < mpKnots2.size(); ++i)
        {
            knot_set2.insert(mpKnots2[i]->Value());
            knots2.push_back(mpKnots2[i]->Value());
        }
        for(std::size_t i = 0; i < mpKnots3.size(); ++i)
        {
            knot_set3.insert(mpKnots3[i]->Value());
            knots3.push_back(mpKnots3[i]->Value());
        }
        
//        std::cout << "knot_set1:";
//        for(std::set<double>::iterator it = knot_set1.begin(); it != knot_set1.end(); ++it)
//            std::cout << " " << *it;
//        std::cout << std::endl;
//        std::cout << "knot_set2:";
//        for(std::set<double>::iterator it = knot_set2.begin(); it != knot_set2.end(); ++it)
//            std::cout << " " << *it;
//        std::cout << std::endl;

        // compute new knots
        std::vector<double> ins_knots1, dknots1(knot_set1.begin(), knot_set1.end());
        std::vector<double> ins_knots2, dknots2(knot_set2.begin(), knot_set2.end());
        std::vector<double> ins_knots3, dknots3(knot_set3.begin(), knot_set3.end());
        if(knots1.size() > 1)
            for(std::size_t i = 0; i < dknots1.size()-1; ++i)
                ins_knots1.push_back(0.5 * (dknots1[i] + dknots1[i+1]));
        if(knots2.size() > 1)
            for(std::size_t i = 0; i < dknots2.size()-1; ++i)
                ins_knots2.push_back(0.5 * (dknots2[i] + dknots2[i+1]));
        if(knots3.size() > 1)
            for(std::size_t i = 0; i < dknots3.size()-1; ++i)
                ins_knots3.push_back(0.5 * (dknots3[i] + dknots3[i+1]));

//        std::cout << "ins_knots1:";
//        for(std::size_t i = 0; i < ins_knots1.size(); ++i)
//            std::cout << " " << ins_knots1[i];
//        std::cout << std::endl;
//        std::cout << "ins_knots2:";
//        for(std::size_t i = 0; i < ins_knots2.size(); ++i)
//            std::cout << " " << ins_knots2[i];
//        std::cout << std::endl;

        // insert the new knot into the knot vector
        for(std::size_t i = 0; i < ins_knots1.size(); ++i)
            for(knot_container_t::iterator it = mpKnots1.begin(); it != mpKnots1.end(); ++it)
                if(ins_knots1[i] < (*it)->Value())
                {
                    knot_t p_knot = knot_t(new Knot(ins_knots1[i]));
                    mpKnots1.insert(it, p_knot);
                    break;
                }
        for(std::size_t i = 0; i < ins_knots2.size(); ++i)
            for(knot_container_t::iterator it = mpKnots2.begin(); it != mpKnots2.end(); ++it)
                if(ins_knots2[i] < (*it)->Value())
                {
                    knot_t p_knot = knot_t(new Knot(ins_knots2[i]));
                    mpKnots2.insert(it, p_knot);
                    break;
                }
        for(std::size_t i = 0; i < ins_knots3.size(); ++i)
            for(knot_container_t::iterator it = mpKnots3.begin(); it != mpKnots3.end(); ++it)
                if(ins_knots3[i] < (*it)->Value())
                {
                    knot_t p_knot = knot_t(new Knot(ins_knots3[i]));
                    mpKnots3.insert(it, p_knot);
                    break;
                }
        
        // update the index of the knots
        for(std::size_t i = 0; i < mpKnots1.size() ; ++i)
            mpKnots1[i]->UpdateIndex(i);
        for(std::size_t i = 0; i < mpKnots2.size() ; ++i)
            mpKnots2[i]->UpdateIndex(i);
        for(std::size_t i = 0; i < mpKnots3.size() ; ++i)
            mpKnots3[i]->UpdateIndex(i);
        
//        std::cout << "mpKnots1:";
//        for(std::size_t i = 0; i < mpKnots1.size() ; ++i)
//            std::cout << " (" << mpKnots1[i]->Index() << "," << mpKnots1[i]->Value() << ")";
//        std::cout << std::endl;
//        std::cout << "mpKnots2:";
//        for(std::size_t i = 0; i < mpKnots2.size() ; ++i)
//            std::cout << " (" << mpKnots2[i]->Index() << "," << mpKnots2[i]->Value() << ")";
//        std::cout << std::endl;

        #ifdef ENABLE_PROFILING
        double time_current = OpenMPUtils::GetCurrentTime();
        std::cout << "Update the knot vectors completed: " << (time_current - time_start) << " s" << std::endl;
        time_start = time_current;
        #endif

        // loop all levels and compute the number of basis functions and cells
        int num_basis_funcs = 0;
        int num_cells = 0;
        for(level_const_iterator it = level_const_begin(); it != level_const_end(); ++it)
        {
            num_basis_funcs += (*it)->NumberOfBasisFunctions();
            num_cells += (*it)->NumberOfCells();
        }

        // create and initialize a new level
        level_t pLevel = level_t(new HnLevel(mpLevels.size() + 1));
        if(mDim == 2)
            pLevel->Initialize(num_cells, num_basis_funcs, mOrder1, mOrder2, mpKnots1, mpKnots2);
        else if(mDim == 3)
            pLevel->Initialize(num_cells, num_basis_funcs, mOrder1, mOrder2, mOrder3, mpKnots1, mpKnots2, mpKnots3);
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
        AddLevel(pLevel);

        #ifdef ENABLE_PROFILING
        time_current = OpenMPUtils::GetCurrentTime();
        std::cout << "Level " << pLevel->Id() << " is initialized: " << (time_current - time_start) << " s" << std::endl;
        time_start = time_current;
        #endif

        // make sure the basis functions on the current level is deactivated
        for(HnLevel::bf_iterator it = pLevel->bf_begin(); it != pLevel->bf_end(); ++it)
            (*it)->Deactivate();

        // compute the refinement coefficient matrix from previous level to new level
        Matrix D;
//        CompressedMatrix D; // will be very slow for deep level
        std::vector<double> new_knots1;
        std::vector<double> new_knots2;
        std::vector<double> new_knots3;

        if(mDim == 2)
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients2D(D,
                                                                     new_knots1,
                                                                     new_knots2,
                                                                     mOrder1,
                                                                     mOrder2,
                                                                     knots1,
                                                                     knots2,
                                                                     ins_knots1,
                                                                     ins_knots2);
        else if(mDim == 3)
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients3D(D,
                                                                     new_knots1,
                                                                     new_knots2,
                                                                     new_knots3,
                                                                     mOrder1,
                                                                     mOrder2,
                                                                     mOrder3,
                                                                     knots1,
                                                                     knots2,
                                                                     knots3,
                                                                     ins_knots1,
                                                                     ins_knots2,
                                                                     ins_knots3);
            
        #ifdef ENABLE_PROFILING
        time_current = OpenMPUtils::GetCurrentTime();
        std::cout << "Compute matrix D completed: " << (time_current - time_start) << " s" << std::endl;
        time_start = time_current;
        #endif

//        KRATOS_WATCH(D)
        // bound check
        level_t pPrevLevel = *(level_end() - 2);
        if(D.size1() != pPrevLevel->NumberOfBasisFunctions())
            KRATOS_THROW_ERROR(std::logic_error, "Something wrong in computing the refinement coefficients. D.size1() must be equal to number of basis functions of previous level", "")
        if(D.size2() != pLevel->NumberOfBasisFunctions())
            KRATOS_THROW_ERROR(std::logic_error, "Something wrong in computing the refinement coefficients. D.size2() must be equal to number of basis functions of this level", "")

        // create linkage from previous level to this level
        double tol = 1.0e-10;
        for(std::size_t i = 0; i < D.size1(); ++i)
            for(std::size_t j = 0; j < D.size2(); ++j)
                if(fabs(D(i, j)) > tol)
                {
                    // get the father basis function
                    HnLevel::bf_iterator it_bf_father = pPrevLevel->bf_begin() + i;

                    // get the child basis function, which the father can be built upon
                    HnLevel::bf_iterator it_bf_child = pLevel->bf_begin() + j;

                    // bind the child to the father
                    (*it_bf_father)->AddChild(*it_bf_child, D(i, j));

                    // connect the cells of the child to the father
                    for(HnBasisFunction::cell_iterator it_cell_father = (*it_bf_father)->cell_begin(); it_cell_father != (*it_bf_father)->cell_end(); ++it_cell_father)
                        for(HnBasisFunction::cell_iterator it_cell_child = (*it_bf_child)->cell_begin(); it_cell_child != (*it_bf_child)->cell_end(); ++it_cell_child)
                        {
//                            std::cout << "children comprises cell " << (*it_cell_child)->Id() << std::endl;
                            if((*it_cell_child)->IsCoverred(*it_cell_father, mDim))
                            {
                                (*it_cell_father)->AddChild(*it_cell_child);
//                                std::cout << "cell " << (*it_cell_father)->Id() << " of bf " << (*it_bf_father)->Id() << " add child " << (*it_cell_child)->Id() << " of bf " << (*it_bf_child)->Id() << ", size = " << (*it_cell_father)->NumberOfChildren() << std::endl;
                            }
                        }
                }

        #ifdef ENABLE_PROFILING
        time_current = OpenMPUtils::GetCurrentTime();
        std::cout << "Build level connection completed: " << (time_current - time_start) << " s" << std::endl;
        time_start = time_current;
        #endif

        std::cout << "Build level " << pLevel->Id() << " completed" << std::endl;
    }

    void HnMesh::RecursiveAddSupportedBasisFunction(cell_t p_cell, bf_t p_bf)
    {
        if(p_cell->IsActive())
        {
            // compute the Bezier extraction of this cell contributing to the basis function
            Vector Crow;
            p_bf->ComputeExtractionOperator(p_cell, Crow, mOrder1, mOrder2);
            p_cell->AddAnchor(p_bf->Id(), p_bf->W(), Crow);
            
            // do a dummy check if this cell contains any active cell. By design, any active cell shall not contain any active children cells.
            for(HnCell::cell_iterator it_cell = p_cell->cell_begin(); it_cell != p_cell->cell_end(); ++it_cell)
                if(it_cell->second->IsActive())
                {
                    std::stringstream ss;
                    ss << "Error detected: active cell shall not contain any active cell. I found cell "
                       << p_cell->Id() << " contains cell " << it_cell->second->Id() << std::endl;
                    ss << "It happened when I searched for supported cells of basis function " << p_bf->Id() << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
        }
        else
        {
            for(HnCell::cell_iterator it_cell = p_cell->cell_begin(); it_cell != p_cell->cell_end(); ++it_cell)
            {
                RecursiveAddSupportedBasisFunction(it_cell->second, p_bf);
            }
        }
    }

    bool HnMesh::IsLinearDependent(bf_t p_bf) const
    {
        // if the bf has no children, which means it can be represented by a linear combination of sub-bfs, then it is not linear dependent
        if(p_bf->NumberOfChildren() == 0)
            return false;
        // otherwise, the bf is linear dependent if all of its children is either active or is represented by a linear combination of sub-bfs
        else
        {
            bool isLinearDependent = true;
            for(HnBasisFunction::bf_const_iterator it_bf = p_bf->bf_const_begin(); it_bf != p_bf->bf_const_end(); ++it_bf)
            {
                if((*it_bf)->IsActive())
                    continue;
                else
                    isLinearDependent = isLinearDependent && IsLinearDependent((*it_bf));
                // if one of the child is not linear dependent, then this bf is also not linear dependent
                if(isLinearDependent == false)
                    return false;
            }
            return true;
        }
    }

    bool HnMesh::IsCoverred(cell_t p_cell) const
    {
        // if the cell has no children, which means it can be represented by a combination of sub-cells, then it is not coverred
        if(p_cell->NumberOfChildren() == 0)
            return false;
        // otherwise, the cell is coverred if all of its children is either active or is represented by a combination of sub-cells
        else
        {
            bool isCoverred = true;
            for(HnCell::cell_const_iterator it_cell = p_cell->cell_const_begin(); it_cell != p_cell->cell_const_end(); ++it_cell)
            {
                if(it_cell->second->IsActive())
                    continue;
                else
                    isCoverred = isCoverred && IsCoverred(it_cell->second);
                // if one of the child is not coverred, then this cell is also not coverred
                if(isCoverred == false)
                    return false;
            }
            return true;
        }
    }

    void HnMesh::RefineWindow(double Xi_min, double Xi_max, double Eta_min, double Eta_max, double Zeta_min, double Zeta_max)
    {
        std::map<unsigned int, std::vector<bf_t> > ToBeRefinedBfs;
        // search and mark all basis functions need to refine on all level (starting from the last level) which support is contained in the refining domain
        double bounding_box_xi_min, bounding_box_xi_max, bounding_box_eta_min, bounding_box_eta_max, bounding_box_zeta_min, bounding_box_zeta_max;
        for(level_iterator it_lvl = level_begin(); it_lvl != level_end(); ++it_lvl)
        {
            unsigned int LevelId = (*it_lvl)->Id();
            // iterate through all active basis function, marking the one need to be refined
            std::vector<bf_t> ToBeRefinedBfsAtLevel;
            for(HnLevel::bf_iterator it_bf = (*it_lvl)->bf_begin(); it_bf != (*it_lvl)->bf_end(); ++it_bf)
            {
                if((*it_bf)->IsActive())
                {
                    // get the bounding box (support domain of the basis function)
                    (*it_bf)->GetBoundingBox(bounding_box_xi_min, bounding_box_xi_max, bounding_box_eta_min, bounding_box_eta_max, bounding_box_zeta_min, bounding_box_zeta_max);

                    // check if the bounding box lie in the refined domain
                    // Remarks: this can be changed by a refinement indicator (i.e from error estimator)
                    if(    bounding_box_xi_min >= Xi_min
                        && bounding_box_xi_max <= Xi_max
                        && bounding_box_eta_min >= Eta_min
                        && bounding_box_eta_max <= Eta_max
                        && bounding_box_zeta_min >= Zeta_min
                        && bounding_box_zeta_max <= Zeta_max )
                    {
                        ToBeRefinedBfsAtLevel.push_back((*it_bf));
                    }
                }
            }
            // store the bfs need to be refined at this level
            if(ToBeRefinedBfsAtLevel.size() != 0)
                ToBeRefinedBfs[LevelId] = ToBeRefinedBfsAtLevel;
        }
        // perform refinement
        Refine(ToBeRefinedBfs);
    }

    void HnMesh::RefineNodes(boost::python::list& pyList)
    {
        // extract the python list to std::set
        std::set<unsigned int> NodeIds;
        typedef boost::python::stl_input_iterator<int> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& node_id, 
                      std::make_pair(iterator_value_type(pyList), // begin
                      iterator_value_type() ) ) // end
            NodeIds.insert(node_id);

//        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
//        {
//            // build a cuboid set from all marked bfs in this level
//            HnCuboid::Pointer p_cuboid_set = HnCuboid::Pointer(new HnCuboid());
//            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
//            {
//                if((*it_bf)->IsActive())
//                {
//                    unsigned int Id = (*it_bf)->Id();
//                    if(std::find(NodeIds.begin(), NodeIds.end(), Id) != NodeIds.end())
//                    {
//                        // get the bounding box (support domain of the basis function)
//                        double bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max;
//                        (*it_bf)->GetBoundingBox(bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max);
//                        p_cuboid_set->AddCuboid(bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max);
//                    }
//                }
//            }
// 
//            // for all the bfs in this level, if it is contained in the cuboid set then it will be added to the refinement list
//            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
//            {
//                if((*it_bf)->IsActive())
//                {
//                    double bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max;
//                    (*it_bf)->GetBoundingBox(bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max);
//                    if(p_cuboid_set->IsInside(bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max))
//                        NodeIds.insert((*it_bf)->Id());
//                }
//            }
//        }

        // now building the list for refinement
        std::map<unsigned int, std::vector<bf_t> > ToBeRefinedBfs;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
        {
            unsigned int LevelId = (*it_lvl)->Id();
            std::vector<bf_t> ToBeRefinedBfsAtLevel;
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
            {
                if((*it_bf)->IsActive())
                {
                    unsigned int Id = (*it_bf)->Id();
                    if(std::find(NodeIds.begin(), NodeIds.end(), Id) != NodeIds.end())
                        ToBeRefinedBfsAtLevel.push_back((*it_bf));
                }
            }
            if(ToBeRefinedBfsAtLevel.size() != 0)
                ToBeRefinedBfs[LevelId] = ToBeRefinedBfsAtLevel;
        }

        // perform the refinement
        this->Refine(ToBeRefinedBfs);
    }

    void HnMesh::Refine(std::map<unsigned int, std::vector<bf_t> >& ToBeRefinedBfs)
    {
        #ifdef DEBUG_REFINE
        for(std::map<unsigned int, std::vector<bf_t> >::iterator it_lvl = ToBeRefinedBfs.begin(); it_lvl != ToBeRefinedBfs.end(); ++it_lvl)
        {
            std::cout << "Bfs will be refined at level " << it_lvl->first << ":";
            for(std::vector<bf_t>::iterator it_bf = it_lvl->second.begin(); it_bf != it_lvl->second.end(); ++it_bf)
                std::cout << " " << (*it_bf)->Id();
            std::cout << std::endl;
        }
        #endif

        while(ToBeRefinedBfs.size() != 0)
        {
            // refine from top level
//            std::set<int> deactivation_list;
            for(std::map<unsigned int, std::vector<bf_t> >::iterator it_lvl = ToBeRefinedBfs.begin(); it_lvl != ToBeRefinedBfs.end(); ++it_lvl)
            {
                #ifdef DEBUG_REFINE
                std::cout << "Refine at level " << it_lvl->first << ":" << std::endl;
                #endif
                std::map<int, HnBasisFunction::Pointer> new_basis_funcs;
                // iterate all bfs in this level
                for(std::vector<bf_t>::iterator it_bf = it_lvl->second.begin(); it_bf != it_lvl->second.end(); ++it_bf)
                {
                    // if the basis function has no children, which means it cannot be refined
                    if((*it_bf)->NumberOfChildren() == 0)
                        continue;

                    #ifdef DEBUG_REFINE
                    std::cout << "bf " << (*it_bf)->Id() << " will be deactivated, bf";
                    #endif

                    // activate the children
                    double father_weight = (*it_bf)->W();
                    double father_X = (*it_bf)->X0();
                    double father_Y = (*it_bf)->Y0();
                    double father_Z = (*it_bf)->Z0();
                    #ifdef DEBUG_REFINE
                    std::set<unsigned int> ActivatedCellIds;
                    std::set<unsigned int> DeactivatedCellIds;
                    #endif
                    for(HnBasisFunction::bf_iterator it_bf_child = (*it_bf)->bf_begin(); it_bf_child != (*it_bf)->bf_end(); ++it_bf_child)
                    {
                        (*it_bf_child)->Activate();

                        #ifdef DEBUG_REFINE
                        std::cout << " " << (*it_bf_child)->Id();
                        #endif

                        double refined_coeff = (*it_bf)->GetRefinedCoefficient((*it_bf_child)->Id());
                        double add_weight = father_weight * refined_coeff;
                        (*it_bf_child)->AddCoordinates(father_X, father_Y, father_Z, add_weight);

                        // activate the cells supporting the child
                        for(HnBasisFunction::cell_iterator it_bf_child_cell = (*it_bf_child)->cell_begin();
                                it_bf_child_cell != (*it_bf_child)->cell_end(); ++it_bf_child_cell)
                        {
//                            bool need_activate = true;
// 
//                            // this child should not be in deactivation_list in order to be activated
//                            need_activate = need_activate && (deactivation_list.find((*it_bf_child_cell)->Id()) == deactivation_list.end());
// 
//                            // none of the children of this cell is active in order to activate the cell, this prevents active cells overlapping (children are active and also father is active)
//                            for(HnCell::cell_iterator it_bf_child_cell_child_cell = (*it_bf_child_cell)->cell_begin();
//                                    it_bf_child_cell_child_cell != (*it_bf_child_cell)->cell_end(); ++it_bf_child_cell_child_cell)
//                                need_activate = need_activate && !(it_bf_child_cell_child_cell->second->IsActive());
// 
//                            if(need_activate)
//                            {
                                (*it_bf_child_cell)->Activate();
                                #ifdef DEBUG_REFINE
                                ActivatedCellIds.insert((*it_bf_child_cell)->Id());
                                #endif
//                            }
                        }

                        // add child to the list of newly activated basis function
                        new_basis_funcs[(*it_bf_child)->Id()] = (*it_bf_child);
                    }

                    #ifdef DEBUG_REFINE
                    std::cout << " is activated" << std::endl;
                    #endif

                    // deactivate the father
                    (*it_bf)->Deactivate();
                    (*it_bf)->SetCoordinates(0.0, 0.0, 0.0, 0.0);
                    
                    // also deactivate the cells supporting the father
                    for(HnBasisFunction::cell_iterator it_bf_cell = (*it_bf)->cell_begin(); it_bf_cell != (*it_bf)->cell_end(); ++it_bf_cell)
                    {
                        (*it_bf_cell)->Deactivate();
                        // put this cell of the father to deactivation list; it will not be activated again anymore
//                        deactivation_list.insert((*it_bf_cell)->Id());
                        #ifdef DEBUG_REFINE
                        DeactivatedCellIds.insert((*it_bf_cell)->Id());
                        #endif
                    }

                    #ifdef DEBUG_REFINE
                    std::cout << "deactivated cells:";
                    for(std::set<unsigned int>::iterator it = DeactivatedCellIds.begin(); it != DeactivatedCellIds.end(); ++it)
                        std::cout << " " << (*it);
                    std::cout << std::endl;

                    std::cout << "activated cells:";
                    for(std::set<unsigned int>::iterator it = ActivatedCellIds.begin(); it != ActivatedCellIds.end(); ++it)
                        std::cout << " " << (*it);
                    std::cout << std::endl;
                    std::cout << "-----------------" << std::endl;
                    #endif
                }

                #ifdef DEBUG_REFINE
                std::cout << "Refine at level " << it_lvl->first << " completed" << std::endl;
                std::cout << "----------------------------" << std::endl;
                #endif
            }
            
            // we search through all level to check for any active bf that are able to be represented by a linear combination of sub-bfs, and we refine it, to ensure a complete linear independence
            ToBeRefinedBfs.clear();
            for(level_iterator it_lvl = level_begin(); it_lvl != level_end(); ++it_lvl)
            {
                unsigned int LevelId = (*it_lvl)->Id();
                std::vector<bf_t> ToBeRefinedBfsAtLevel;
                for(HnLevel::bf_iterator it_bf = (*it_lvl)->bf_begin(); it_bf != (*it_lvl)->bf_end(); ++it_bf)
                    if((*it_bf)->IsActive())
                        if(IsLinearDependent((*it_bf)))
                            ToBeRefinedBfsAtLevel.push_back((*it_bf));
                if(ToBeRefinedBfsAtLevel.size() != 0)
                    ToBeRefinedBfs[LevelId] = ToBeRefinedBfsAtLevel;
                #ifdef DEBUG_REFINE
                std::cout << "Level " << LevelId << " needs to refine additional " << ToBeRefinedBfsAtLevel.size() << " bfs:";
                for(unsigned int i = 0; i < ToBeRefinedBfsAtLevel.size(); ++i)
                    std::cout << " " << ToBeRefinedBfsAtLevel[i]->Id();
                std::cout << std::endl;
                #endif
            }
            #ifdef DEBUG_REFINE
            std::cout << "Refinement cycle completed" << std::endl;
            std::cout << "----------------------------------" << std::endl;
            #endif
        }

        // we seach all the cells, to see if it is coverred by a combination of sub-cells, if it is, then deactivate it
        for(level_iterator it_lvl = level_begin(); it_lvl != level_end(); ++it_lvl)
            for(HnLevel::cell_iterator it_cell = (*it_lvl)->cell_begin(); it_cell != (*it_lvl)->cell_end(); ++it_cell)
                if(it_cell->second->IsActive())
                    if(IsCoverred(it_cell->second))
                        it_cell->second->Deactivate();

        // clear all internal data on supported basis functions of the cells
        for(level_iterator it_lvl = level_begin(); it_lvl != level_end(); ++it_lvl)
            for(HnLevel::cell_iterator it_cell = (*it_lvl)->cell_begin(); it_cell != (*it_lvl)->cell_end(); ++it_cell)
                it_cell->second->Reset();

        // for all the cells, check all the supported basis functions to this cell and add to list
        for(level_iterator it_lvl = level_begin(); it_lvl != level_end(); ++it_lvl)
            for(HnLevel::bf_iterator it_bf = (*it_lvl)->bf_begin(); it_bf != (*it_lvl)->bf_end(); ++it_bf)
            {
                if((*it_bf)->IsActive())
                {
                    for(HnBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                    {
                        // add the supported basis functions and compute the Bezier extraction operator over cells
                        RecursiveAddSupportedBasisFunction(*it_cell, *it_bf);
                    }
                }
            }

        std::cout << "Refine completed" << std::endl;
    }

    void HnMesh::ExportCellTopology(std::string fn, bool cell_numbering) const
    {
        std::ofstream outfile(fn.c_str());
        outfile << "%% Cell topology generated from hierarchical NURBS mesh, (c) Hoang Giang Bui, 2015\n";
        outfile << "clc\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";

        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::cell_const_iterator it_cell = (*it_lvl)->cell_const_begin(); it_cell != (*it_lvl)->cell_const_end(); ++it_cell)
            {
                HnLevel::cell_t p_cell = it_cell->second;
                if(p_cell->IsActive())
//                    outfile << *p_cell << std::endl;
                    if(mDim == 2)
                    {
                        outfile << "line([" << p_cell->LeftValue() << " " << p_cell->RightValue() << "],[" << p_cell->DownValue() << " " << p_cell->DownValue() << "]);\n";
                        outfile << "line([" << p_cell->LeftValue() << " " << p_cell->RightValue() << "],[" << p_cell->UpValue() << " " << p_cell->UpValue() << "]);\n";
                        outfile << "line([" << p_cell->LeftValue() << " " << p_cell->LeftValue() << "],[" << p_cell->DownValue() << " " << p_cell->UpValue() << "]);\n";
                        outfile << "line([" << p_cell->RightValue() << " " << p_cell->RightValue() << "],[" << p_cell->DownValue() << " " << p_cell->UpValue() << "]);\n";
                        if(cell_numbering)
                            outfile << "text(" << 0.5*(p_cell->LeftValue() + p_cell->RightValue()) << "," << 0.5*(p_cell->DownValue() + p_cell->UpValue()) << ",'" << p_cell->Id() << "');\n";
                    }
                    else if(mDim == 3)
                    {
                        //TODO
                        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not yet implemented in 3D")
                    }
            }

        outfile.close();
        
        std::cout << "Export cell topology to " << fn << " completed" << std::endl;
    }

    void HnMesh::ExportBasisFunctionTopology(std::string fn, bool bf_numbering, boost::python::list& pyList) const
    {
        // extract the python list to std::set
        std::set<unsigned int> BfIds;
        typedef boost::python::stl_input_iterator<int> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& bf_id, 
                      std::make_pair(iterator_value_type(pyList), // begin
                      iterator_value_type() ) ) // end
            BfIds.insert(bf_id);

        std::ofstream outfile(fn.c_str());
        outfile << "%% Basis function topology generated from hierarchical NURBS mesh, (c) Hoang Giang Bui, 2015\n";
        outfile << "clc\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";

        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
        {
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
            {
                if(std::find(BfIds.begin(), BfIds.end(), (*it_bf)->Id()) != BfIds.end())
                {
                    // get the bounding box of the basis function
                    double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
                    (*it_bf)->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);

                    if(mDim == 2)
                    {
                        outfile << "line([" << Xmin << " " << Xmax << "],[" << Ymin << " " << Ymin << "]);\n";
                        outfile << "line([" << Xmin << " " << Xmax << "],[" << Ymax << " " << Ymax << "]);\n";
                        outfile << "line([" << Xmin << " " << Xmin << "],[" << Ymin << " " << Ymax << "]);\n";
                        outfile << "line([" << Xmax << " " << Xmax << "],[" << Ymin << " " << Ymax << "]);\n";
                        if(bf_numbering)
                            outfile << "text(" << 0.5*(Xmin + Xmax) << "," << 0.5*(Ymin + Ymax) << ",'" << (*it_lvl)->Id() << "," << (*it_bf)->Id() << "');\n";
                    }
                    else if(mDim == 3)
                    {
                        //TODO
                        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not yet implemented in 3D")
                    }
                }
            }
        }

        outfile.close();

        std::cout << "Export basis function topology to " << fn << " completed" << std::endl;
    }

    void HnMesh::GenerateCellGeology(std::vector<unsigned int>& point_list,
                                     std::map<unsigned int, double>& X_list,
                                     std::map<unsigned int, double>& Y_list,
                                     std::map<unsigned int, double>& Z_list,
                                     std::map<unsigned int, double>& xi_list,
                                     std::map<unsigned int, double>& eta_list,
                                     std::map<unsigned int, double>& zeta_list,
                                     std::map<unsigned int, unsigned int>& cell_list,
                                     std::set<unordered_pair<unsigned int> >& Lines,
                                     std::set<unordered_tuple<unsigned int> >& Faces,
                                     std::map<unsigned int, std::vector<std::vector<unsigned int> > >& Connectivities) const
    {
        // get the list of active cells in the hierarchical mesh
        std::vector<cell_t> active_cells;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::cell_const_iterator it_cell = (*it_lvl)->cell_const_begin(); it_cell != (*it_lvl)->cell_const_end(); ++it_cell)
                if(it_cell->second->IsActive())
                    active_cells.push_back(it_cell->second);

        // get the list of active basis functions in the hierarchical mesh
        std::map<unsigned int, bf_t> active_bfs;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
                if((*it_bf)->IsActive())
                    active_bfs[(*it_bf)->Id()] = (*it_bf);

        // extract all the cell's vertices and lines
        AutoCollapseSpatialBinning Binning(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0e-6);
        std::map<unsigned int, cell_t> MapVertexToCell;
        std::map<unsigned int, std::vector<unsigned int> > MapCellToNodes;
        for(std::vector<cell_t>::const_iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            if(mDim == 2)
            {
                unsigned int V1 = Binning.AddNode((*it)->LeftValue(), (*it)->DownValue(), 0.0);
                unsigned int V2 = Binning.AddNode((*it)->RightValue(), (*it)->DownValue(), 0.0);
                unsigned int V3 = Binning.AddNode((*it)->LeftValue(), (*it)->UpValue(), 0.0);
                unsigned int V4 = Binning.AddNode((*it)->RightValue(), (*it)->UpValue(), 0.0);

                MapVertexToCell[V1] = (*it);
                MapVertexToCell[V2] = (*it);
                MapVertexToCell[V3] = (*it);
                MapVertexToCell[V4] = (*it);
                MapCellToNodes[(*it)->Id()].push_back(V1);
                MapCellToNodes[(*it)->Id()].push_back(V2);
                MapCellToNodes[(*it)->Id()].push_back(V4);
                MapCellToNodes[(*it)->Id()].push_back(V3);

                Lines.insert(unordered_pair<unsigned int>(V1, V2));
                Lines.insert(unordered_pair<unsigned int>(V2, V4));
                Lines.insert(unordered_pair<unsigned int>(V4, V3));
                Lines.insert(unordered_pair<unsigned int>(V3, V1));
            }
        }

        // compute the global coordinates of the vertex
        for(unsigned int i = 0; i < Binning.NumberOfNodes(); ++i)
        {
            unsigned int Id = i + 1;
            cell_t p_cell = MapVertexToCell[Id];
            
            if(mDim == 2)
            {
                // compute the local coordinates of the vertex
                double xi = Binning.GetX(Id);
                double eta = Binning.GetY(Id);
                double local_xi = (xi - p_cell->LeftValue()) / (p_cell->RightValue() - p_cell->LeftValue());
                double local_eta = (eta - p_cell->DownValue()) / (p_cell->UpValue() - p_cell->DownValue());

                // compute the Bernstein basis function on the local coordinates
                Vector B((mOrder1 + 1) * (mOrder2 + 1));
                for(unsigned int k = 0; k < mOrder1 + 1; ++k)
                    for(unsigned int l = 0; l < mOrder2 + 1; ++l)
                    {
                        unsigned int num = k * (mOrder2 + 1) + l;
                        double B1 = BezierUtils::bernstein2(k, mOrder1, local_xi);
                        double B2 = BezierUtils::bernstein2(l, mOrder2, local_eta);
                        B(num) = B1 * B2;
                    }

                // get the extraction operator on this cell
                Matrix C = p_cell->GetExtractionOperator();

                // compute the B-splines shape function values at the local coordinates
                Vector N(p_cell->NumberOfAnchors());
                noalias(N) = prod(C, B);

                // compute the NURBS shape function at the local coordinates
                Vector W;
                p_cell->GetAnchorWeights(W);
                double Denom = inner_prod(W, N);
                Vector R(p_cell->NumberOfAnchors());
                for(std::size_t k = 0; k < p_cell->NumberOfAnchors(); ++k)
                    R(k) = W(k) * N(k) / Denom;

                // get the list of supported basis function of this cell
                const std::vector<unsigned int>& bfs_id = p_cell->GetSupportedAnchors();

                // compute the glocal coordinates at the local coordinates
                double X = 0.0, Y = 0.0, Z = 0.0;
                for(std::size_t k = 0; k < p_cell->NumberOfAnchors(); ++k)
                {
                    X += R(k) * active_bfs[bfs_id[k]]->X0();
                    Y += R(k) * active_bfs[bfs_id[k]]->Y0();
                }

                // add point to list
                X_list[Id] = X;
                Y_list[Id] = Y;
                Z_list[Id] = 0.0;
                point_list.push_back(Id);
                xi_list[Id] = local_xi;
                eta_list[Id] = local_eta;
                zeta_list[Id] = 0;
                cell_list[Id] = p_cell->Id();
            }
        }

        // check if any line/face contain any vertex and delete that line
        for(std::vector<cell_t>::const_iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            if(mDim == 2)
            {
                // get the lines of the cell
                unsigned int e1, e2;
                std::vector<unsigned int> MiddleNodes;
                for(std::size_t i = 0; i < MapCellToNodes[(*it)->Id()].size(); ++i)
                {
                    e1 = MapCellToNodes[(*it)->Id()][i];
                    if(i == MapCellToNodes[(*it)->Id()].size() - 1)
                        e2 = MapCellToNodes[(*it)->Id()][0];
                    else
                        e2 = MapCellToNodes[(*it)->Id()][i+1];
                    if(e1 > e2)
                    {
                        unsigned int tmp = e1;
                        e1 = e2;
                        e2 = tmp;
                    }

                    double xi1 = Binning.GetX(e1);
                    double eta1 = Binning.GetY(e1);
                    double xi2 = Binning.GetX(e2);
                    double eta2 = Binning.GetY(e2);
                    if(xi1 > xi2)
                    {
                        double tmp = xi1;
                        xi1 = xi2;
                        xi2 = tmp;
                    }
                    if(eta1 > eta2)
                    {
                        double tmp = eta1;
                        eta1 = eta2;
                        eta2 = tmp;
                    }
//                    KRATOS_WATCH(xi1)
//                    KRATOS_WATCH(xi2)
//                    KRATOS_WATCH(eta1)
//                    KRATOS_WATCH(eta2)
                    // check if the line contain any vertex
                    for(unsigned int i = 0; i < Binning.NumberOfNodes(); ++i)
                    {
                        unsigned int Id = i + 1;
                        double xi = Binning.GetX(Id);
                        double eta = Binning.GetY(Id);
                        if((eta==eta1 && eta==eta2 && xi1<xi && xi<xi2)
                        || (xi==xi1 && xi==xi2 && eta1<eta && eta<eta2))
                        {
//                            KRATOS_WATCH(Id)
                            MiddleNodes.push_back(Id);

                            // find the line which is e1,e2 and delete
                            std::set<unordered_pair<unsigned int> >::iterator it_to_be_deleted = Lines.end();
                            for(std::set<unordered_pair<unsigned int> >::iterator it = Lines.begin(); it != Lines.end(); ++it)
                                if(it->first() == e1 && it->second() == e2)
                                {
                                    it_to_be_deleted = it;
                                    break;
                                }
                            if(it_to_be_deleted != Lines.end())
                                Lines.erase(it_to_be_deleted);
                        }
                    }
                }

                if(MiddleNodes.size() == 0)
                {
                    // the cell has 4 vertices, now generate the connectivities
                    std::vector<unsigned int> Conn(MapCellToNodes[(*it)->Id()].begin(), MapCellToNodes[(*it)->Id()].end());
                    Connectivities[(*it)->Id()].push_back(Conn);
                }
                else
                {
                    // insert the middle vertices to the cell
                    MapCellToNodes[(*it)->Id()].insert( MapCellToNodes[(*it)->Id()].end(),
                                                        MiddleNodes.begin(),
                                                        MiddleNodes.end() );

                    // prepare the list of points
                    std::vector<double> XYlist;
                    for(std::size_t i = 0; i < MapCellToNodes[(*it)->Id()].size(); ++i)
                    {
                        XYlist.push_back(Binning.GetX(MapCellToNodes[(*it)->Id()][i]));
                        XYlist.push_back(Binning.GetY(MapCellToNodes[(*it)->Id()][i]));
                    }

                    // generate the triangulation
                    TriangulationUtils util;
                    std::vector<std::vector<unsigned int> > tris;
                    util.ComputeDelaunayTriangulation(XYlist, tris);

                    // map the triangle connectivities to node Id
                    std::vector<std::vector<unsigned int> > Conns(tris.size());
                    for(std::size_t i = 0; i < tris.size(); ++i)
                        for(std::size_t j = 0; j < tris[i].size(); ++j)
                            Conns[i].push_back(MapCellToNodes[(*it)->Id()][tris[i][j]]);

                    Connectivities[(*it)->Id()] = Conns;
                }
            }
            else if(mDim == 3)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The algorithm in 3D is not yet implemented", "")
            }
        }
    }

    void HnMesh::ExportCellGeology(std::string fn) const
    {
        // generate cell geology
        std::vector<unsigned int> point_list;
        std::map<unsigned int, double> X_list;
        std::map<unsigned int, double> Y_list;
        std::map<unsigned int, double> Z_list;
        std::map<unsigned int, double> xi_list;
        std::map<unsigned int, double> eta_list;
        std::map<unsigned int, double> zeta_list;
        std::map<unsigned int, unsigned int> cell_list;
        std::set<unordered_pair<unsigned int> > Lines;
        std::set<unordered_tuple<unsigned int> > Faces;
        std::map<unsigned int, std::vector<std::vector<unsigned int> > > Connectivities;
        GenerateCellGeology(point_list, X_list, Y_list, Z_list, xi_list, eta_list, zeta_list, cell_list, Lines, Faces, Connectivities);

        // create file handler
        std::ofstream outfile(fn.c_str());
        outfile << "%% Hierarchical NURBS mesh cell geology, (c) Hoang Giang Bui, 2015\n";
        outfile << "clc\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";
        if(mDim == 2)
            for(std::set<unordered_pair<unsigned int> >::iterator it = Lines.begin(); it != Lines.end(); ++it)
                outfile << "line([" << X_list[it->first()] << " " << X_list[it->second()] << "],[" << Y_list[it->first()] << " " << Y_list[it->second()] << "]);\n";

        outfile.close();

        std::cout << "Export cell geology to " << fn << " completed" << std::endl;
    }

    void HnMesh::ExportCellGeologyAsPostMDPA(std::string fn) const
    {
        // generate cell geology
        std::vector<unsigned int> point_list;
        std::map<unsigned int, double> X_list;
        std::map<unsigned int, double> Y_list;
        std::map<unsigned int, double> Z_list;
        std::map<unsigned int, double> xi_list;
        std::map<unsigned int, double> eta_list;
        std::map<unsigned int, double> zeta_list;
        std::map<unsigned int, unsigned int> cell_list;
        std::set<unordered_pair<unsigned int> > Lines;
        std::set<unordered_tuple<unsigned int> > Faces;
        std::map<unsigned int, std::vector<std::vector<unsigned int> > > connectivities_list;
        GenerateCellGeology(point_list, X_list, Y_list, Z_list, xi_list, eta_list, zeta_list, cell_list, Lines, Faces, connectivities_list);
        
        // export to post MDPA
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for post processing of Bezier-based discretization\n";
        outfile << "//(c) 2015 Hoang Giang Bui, Ruhr-University Bochum\n";

        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        outfile << "//This file is created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900)
                << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl << std::endl;

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " " << X_list[point_list[i]]
                                     << " " << Y_list[point_list[i]]
                                     << " " << Z_list[point_list[i]]
                                     << "\n";
        outfile << "End Nodes\n\n";

        // write elements
        if(mDim == 2)
        {
            unsigned int ElementId = 0;

            // write triangle elements
            outfile << "Begin Elements KinematicLinear2D3N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 3)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";

            // write quadrilateral elements
            outfile << "Begin Elements KinematicLinear2D4N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 4)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";
        }
        else if(mDim == 3)
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not yet to be implemented in 3D")
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")
        
        // write nodal data
        outfile << "Begin NodalData LOCAL_COORDINATES\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 [3] (" << xi_list[point_list[i]]
                                     << "," << eta_list[point_list[i]]
                                     << "," << zeta_list[point_list[i]] << ")\n";
        outfile << "End NodalData\n\n";

        outfile << "Begin NodalData PARENT_ELEMENT_ID\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 " << cell_list[point_list[i]] << "\n";
        outfile << "End NodalData\n\n";
        
        std::cout << "Export post MDPA to " << fn << " completed" << std::endl;
    }

    void HnMesh::ExportMatlab(std::string fn) const
    {
        std::ofstream outfile(fn.c_str());
        outfile << "%% Hierarchical NURBS mesh information, (c) Hoang Giang Bui, 2015\n";
        outfile << "clc\n";
        outfile << "clear\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n\n";

        // export the basis function information
        int cnt = 0;
        std::map<int, int> bfMap; // map from bf Id to the local id
        std::vector<double> LocalKnots1, LocalKnots2, LocalKnots3;
        double min_xi = static_cast<double>(INT_MAX);
        double max_xi = -min_xi;
        double min_eta = min_xi;
        double max_eta = -min_eta;
        double min_zeta = min_xi;
        double max_zeta = -min_zeta;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_begin(); it_bf != (*it_lvl)->bf_end(); ++it_bf)
                if((*it_bf)->IsActive())
                {
                    (*it_bf)->GetLocalKnots(1, LocalKnots1);
                    (*it_bf)->GetLocalKnots(2, LocalKnots2);
                    if(mDim == 3)
                        (*it_bf)->GetLocalKnots(3, LocalKnots3);

                    double min_xi_bf = *std::min_element(LocalKnots1.begin(), LocalKnots1.end());
                    double max_xi_bf = *std::max_element(LocalKnots1.begin(), LocalKnots1.end());
                    if(min_xi_bf < min_xi)
                        min_xi = min_xi_bf;
                    if(max_xi_bf > max_xi)
                        max_xi = max_xi_bf;

                    double min_eta_bf = *std::min_element(LocalKnots2.begin(), LocalKnots2.end());
                    double max_eta_bf = *std::max_element(LocalKnots2.begin(), LocalKnots2.end());
                    if(min_eta_bf < min_eta)
                        min_eta = min_eta_bf;
                    if(max_eta_bf > max_eta)
                        max_eta = max_eta_bf;

                    if(mDim == 3)
                    {
                        double min_zeta_bf = *std::min_element(LocalKnots3.begin(), LocalKnots3.end());
                        double max_zeta_bf = *std::max_element(LocalKnots3.begin(), LocalKnots3.end());
                        if(min_zeta_bf < min_zeta)
                            min_zeta = min_zeta_bf;
                        if(max_zeta_bf > max_zeta)
                            max_zeta = max_zeta_bf;
                    }

                    ++cnt;
                    bfMap[(*it_bf)->Id()] = cnt;

                    outfile << "Xi{" << cnt << "} = [";
                    for(std::size_t i = 0; i < LocalKnots1.size(); ++i)
                        outfile << " " << LocalKnots1[i];
                    outfile << "];\n";

                    outfile << "Eta{" << cnt << "} = [";
                    for(std::size_t i = 0; i < LocalKnots2.size(); ++i)
                        outfile << " " << LocalKnots2[i];
                    outfile << "];\n";

                    if(mDim == 3)
                    {
                        outfile << "Zeta{" << cnt << "} = [";
                        for(std::size_t i = 0; i < LocalKnots3.size(); ++i)
                            outfile << " " << LocalKnots3[i];
                        outfile << "];\n";
                    }

                    outfile << "P(" << cnt << ",:) = [" << (*it_bf)->X0() << " " << (*it_bf)->Y0() << " " << (*it_bf)->Z0() << "];\n";
                    outfile << "W(" << cnt << ") = " << (*it_bf)->W() << ";\n";
                    outfile << "Id(" << cnt << ") = " << (*it_bf)->Id() << ";\n";
//                    for(HnBasisFunction::cell_const_iterator it_cell = (*it_bf)->cell_const_begin(); it_cell != (*it_bf)->cell_const_end(); ++it_cell)
//                        outfile << "%" << *(*it_cell) << std::endl;
                    outfile << std::endl;
                }

        // export the cell information
        // get the list of active cells in the hierarchical mesh
        std::vector<cell_t> active_cells;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::cell_const_iterator it_cell = (*it_lvl)->cell_const_begin(); it_cell != (*it_lvl)->cell_const_end(); ++it_cell)
                if(it_cell->second->IsActive())
                    active_cells.push_back(it_cell->second);

        cnt = 0;
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            ++cnt;

            // write the boundary of the cell
            outfile << "S{" << cnt << "} = [" << (*it)->LeftValue() << " " << (*it)->RightValue() << ";";
            outfile << (*it)->DownValue() << " " << (*it)->UpValue() << "];\n";

            // write the extraction operator
            Matrix C = (*it)->GetExtractionOperator();
            outfile << "C{" << cnt << "} = [";
            for(std::size_t i = 0; i < C.size1(); ++i)
            {
                for(std::size_t j = 0;  j < C.size2(); ++ j)
                    outfile << " " << C(i, j);
                outfile << ";";
            }
            outfile << "];\n";
            
            // write the supported basis functions
            const std::vector<unsigned int>& bfs = (*it)->GetSupportedAnchors();
            outfile << "N{" << cnt << "} = [";
            for(std::size_t i = 0; i < bfs.size(); ++i)
                outfile << " " << bfs[i];
            outfile << "];\n" << std::endl;
        }

        outfile << std::endl;

        outfile.close();
        std::cout << "Export mesh information to " << fn << " completed" << std::endl;
    }

    void HnMesh::ExportMDPA(std::string fn) const
    {
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for hierarchical NURBS\n";
        outfile << "//(c) 2015 Hoang Giang Bui, Ruhr-University Bochum\n";
        
        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        outfile << "//This file is created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900)
                << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl << std::endl;

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        int BfId = 0;
        #ifdef MDPA_NODE_RENUMBERING
        std::map<int, int> OldToNewBfId; // map from old basis function id to new basis function id
        #endif
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_begin(); it_bf != (*it_lvl)->bf_end(); ++it_bf)
                if((*it_bf)->IsActive())
                {
                    #ifdef MDPA_NODE_RENUMBERING
                    outfile << ++BfId << "\t" << (*it_bf)->X0() << "\t" << (*it_bf)->Y0() << "\t" << (*it_bf)->Z0() << std::endl; // assign new basis function id
                    OldToNewBfId[(*it_bf)->Id()] = BfId;
                    #else
                    outfile << (*it_bf)->Id() << "\t" << (*it_bf)->X0() << "\t" << (*it_bf)->Y0() << "\t" << (*it_bf)->Z0() << std::endl; // assign new basis function id
                    #endif
                }
        outfile << "End Nodes\n\n";

        // get the list of active cells in the hierarchical mesh
        std::vector<cell_t> active_cells;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::cell_const_iterator it_cell = (*it_lvl)->cell_const_begin(); it_cell != (*it_lvl)->cell_const_end(); ++it_cell)
                if(it_cell->second->IsActive())
                    active_cells.push_back(it_cell->second);

        // write elements
        if(mDim == 2)
            outfile << "Begin Elements KinematicLinearGeo2dBezier\n";
        else if(mDim == 3)
            outfile << "Begin Elements KinematicLinearGeo3dBezier\n";
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")
        #ifdef MDPA_CELL_RENUMBERING
        int ElemId = 0;
        std::map<int, int> OldToNewElemId; // map from old element id to new element id
        #endif
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            const std::vector<unsigned int>& supported_bfs = (*it)->GetSupportedAnchors();
            #ifdef MDPA_CELL_RENUMBERING
            outfile << ++ElemId << " 1"; // assign new element Id
            OldToNewElemId[(*it)->Id()] = ElemId;
            #else
            outfile << (*it)->Id() << " 1"; // assign new element Id
            #endif
            for(std::size_t i = 0; i < supported_bfs.size(); ++i)
                #ifdef MDPA_NODE_RENUMBERING
                outfile << " " << OldToNewBfId[supported_bfs[i]];
                #else
                outfile << " " << supported_bfs[i];
                #endif
            outfile << std::endl;
        }
        outfile << "End Elements\n\n";

        // write weights
        outfile << "Begin ElementalData NURBS_WEIGHT\n";
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            const std::vector<double>& weights = (*it)->GetAnchorWeights();
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " [" << weights.size() << "] (";
            #else
            outfile << (*it)->Id() << " [" << weights.size() << "] (";
            #endif
            for(std::size_t i = 0; i < weights.size() - 1; ++i)
                outfile << weights[i] << ",";
            outfile << weights.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // read extraction operator from element in the form of triplet CSR
        // this requires that the Id of the cell must be unique
        std::vector<int> rowPtr;
        std::vector<int> colInd;
        std::vector<double> values;
        std::map<int, std::vector<int> > rowPtrMap;
        std::map<int, std::vector<int> > colIndMap;
        std::map<int, std::vector<double> > valuesMap;
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
        {
            rowPtr.clear();
            colInd.clear();
            values.clear();
            (*it)->GetExtractionOperator(rowPtr, colInd, values);
            rowPtrMap[(*it)->Id()] = rowPtr;
            colIndMap[(*it)->Id()] = colInd;
            valuesMap[(*it)->Id()] = values;
        }

        // write extraction operator
        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_ROWPTR\n";
        for(std::map<int, std::vector<int> >::iterator it = rowPtrMap.begin(); it != rowPtrMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_COLIND\n";
        for(std::map<int, std::vector<int> >::iterator it = colIndMap.begin(); it != colIndMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_VALUES\n";
        for(std::map<int, std::vector<double> >::iterator it = valuesMap.begin(); it != valuesMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // write the degree
        outfile << "Begin ElementalData NURBS_DEGREE_1\n";
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder1 << std::endl;
            #else
            outfile << (*it)->Id() << " " << mOrder1 << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData NURBS_DEGREE_2\n";
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder2 << std::endl;
            #else
            outfile << (*it)->Id() << " " << mOrder2 << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        if(mDim == 3)
        {
            outfile << "Begin ElementalData NURBS_DEGREE_3\n";
            for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
                #ifdef MDPA_CELL_RENUMBERING
                outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder3 << std::endl;
                #else
                outfile << (*it)->Id() << " " << mOrder3 << std::endl;
                #endif
            outfile << "End ElementalData\n\n";
        }

        // write the division
        outfile << "Begin ElementalData NUM_DIVISION_1\n";
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
            #else
            outfile << (*it)->Id() << " 1" << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData NUM_DIVISION_2\n";
        for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
            #else
            outfile << (*it)->Id() << " 1" << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        if(mDim == 3)
        {
            outfile << "Begin ElementalData NUM_DIVISION_3\n";
            for(std::vector<cell_t>::iterator it = active_cells.begin(); it != active_cells.end(); ++it)
                #ifdef MDPA_CELL_RENUMBERING
                outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
                #else
                outfile << (*it)->Id() << " 1" << std::endl;
                #endif
            outfile << "End ElementalData\n\n";
        }

        outfile.close();
        std::cout << "Export MDPA to " << fn << " completed" << std::endl;
    }

    void HnMesh::ExportPostMDPA(std::string fn, int NumDivision1, int NumDivision2, int NumDivision3) const
    {
        // create a vector of sampling knots in each direction
        std::vector<double> SamplingKnots1;
        std::vector<double> SamplingKnots2;
        std::vector<double> SamplingKnots3;

        int min_xi = (*(mpKnots1.begin()))->Value();
        int max_xi = (*(mpKnots1.end() - 1))->Value();
        for(std::size_t i = 0; i < NumDivision1 + 1; ++i)
            SamplingKnots1.push_back(min_xi + (double) i / NumDivision1 * (max_xi - min_xi));

        int min_eta = (*(mpKnots2.begin()))->Value();
        int max_eta = (*(mpKnots2.end() - 1))->Value();
        for(std::size_t i = 0; i < NumDivision2 + 1; ++i)
            SamplingKnots2.push_back(min_eta + (double) i / NumDivision2 * (max_eta - min_eta));

        if(mDim == 3)
        {
            int min_zeta = (*(mpKnots3.begin()))->Value();
            int max_zeta = (*(mpKnots3.end() - 1))->Value();
            for(std::size_t i = 0; i < NumDivision3 + 1; ++i)
                SamplingKnots3.push_back(min_zeta + (double) i / NumDivision3 * (max_zeta - min_zeta));
        }

        // get the list of active cells in the hierarchical mesh
        std::vector<cell_t> active_cells;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::cell_const_iterator it_cell = (*it_lvl)->cell_const_begin(); it_cell != (*it_lvl)->cell_const_end(); ++it_cell)
                if(it_cell->second->IsActive())
                    active_cells.push_back(it_cell->second);

        // get the list of active basis functions in the hierarchical mesh
        std::map<int, bf_t> active_bfs;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
                if((*it_bf)->IsActive())
                    active_bfs[(*it_bf)->Id()] = (*it_bf);

        // generate a list of point based on sampling points
        std::vector<int> point_list;
        std::map<int, double> X_list;
        std::map<int, double> Y_list;
        std::map<int, double> Z_list;
        std::map<int, double> xi_list;
        std::map<int, double> eta_list;
        std::map<int, double> zeta_list;
        std::map<int, int> parent_element_list;
        std::vector<std::vector<unsigned int> > connectivities_list;
        if(mDim == 2)
        {
            int cnt = 0;
            Vector W;
            for(std::size_t i = 0; i < SamplingKnots1.size(); ++i)
            {
                for(std::size_t j = 0; j < SamplingKnots2.size(); ++j)
                {
                    // check if this point lie in which cell
                    bool found = false;
                    for(std::vector<cell_t>::iterator it_cell = active_cells.begin(); it_cell != active_cells.end(); ++it_cell)
                    {
                        if((*it_cell)->IsCoverage(SamplingKnots1[i], SamplingKnots2[j]))
                        {
                            found = true;

                            // compute the local coordinates of the sampling point
                            double local_xi = (SamplingKnots1[i] - (*it_cell)->LeftValue()) / ((*it_cell)->RightValue() - (*it_cell)->LeftValue());
                            double local_eta = (SamplingKnots2[j] - (*it_cell)->DownValue()) / ((*it_cell)->UpValue() - (*it_cell)->DownValue());
//                            KRATOS_WATCH(local_xi)
//                            KRATOS_WATCH(local_eta)

                            // compute the Bernstein basis function on the local coordinates
                            Vector B((mOrder1 + 1) * (mOrder2 + 1));
                            for(unsigned int k = 0; k < mOrder1 + 1; ++k)
                                for(unsigned int l = 0; l < mOrder2 + 1; ++l)
                                {
                                    unsigned int num = k * (mOrder2 + 1) + l;
                                    double B1 = BezierUtils::bernstein2(k, mOrder1, local_xi);
                                    double B2 = BezierUtils::bernstein2(l, mOrder2, local_eta);
                                    B(num) = B1 * B2;
                                }

                            // get the extraction operator on this cell
                            Matrix C = (*it_cell)->GetExtractionOperator();

                            // compute the B-splines shape function values at the local coordinates
                            Vector N((*it_cell)->NumberOfAnchors());
                            noalias(N) = prod(C, B);

                            // compute the NURBS shape function at the local coordinates
                            (*it_cell)->GetAnchorWeights(W);
                            double Denom = inner_prod(W, N);
                            Vector R((*it_cell)->NumberOfAnchors());
                            for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                R(k) = W(k) * N(k) / Denom;

                            // get the list of supported basis function of this cell
                            const std::vector<unsigned int>& bfs_id = (*it_cell)->GetSupportedAnchors();

                            // compute the glocal coordinates at the local coordinates
                            double X = 0.0, Y = 0.0, Z = 0.0;
                            for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                            {
                                X += R(k) * active_bfs[bfs_id[k]]->X0();
                                Y += R(k) * active_bfs[bfs_id[k]]->Y0();
                                Z += R(k) * active_bfs[bfs_id[k]]->Z0();
                            }

                            // add point to list
                            point_list.push_back(++cnt);
                            X_list[cnt] = X;
                            Y_list[cnt] = Y;
                            Z_list[cnt] = Z;

                            xi_list[cnt] = local_xi;
                            eta_list[cnt] = local_eta;
                            zeta_list[cnt] = 0.0;

                            parent_element_list[cnt] = (*it_cell)->Id();
                            break;
                        }
                    }
                    if(found == false)
                        KRATOS_THROW_ERROR(std::logic_error, "The sampling points is not detected in the hierarchical mesh", "")
                }
            }

            // generate the element connectivities
            std::vector<unsigned int> connectivities(4);
            for(std::size_t i = 0; i < NumDivision1; ++i)
            {
                for(std::size_t j = 0; j < NumDivision2; ++j)
                {
//                    connectivities[0] = i * (NumDivision2 + 1) + j + 1;
//                    connectivities[1] = connectivities[0] + 1;
//                    connectivities[3] = connectivities[0] + NumDivision2 + 1;
//                    connectivities[2] = connectivities[3] + 1;
//                    connectivities_list.push_back(connectivities);

                    unsigned int Node1 = i * (NumDivision2 + 1) + j + 1;
                    unsigned int Node2 = i * (NumDivision2 + 1) + j + 2;
                    unsigned int Node3 = (i + 1) * (NumDivision2 + 1) + j + 1;
                    unsigned int Node4 = (i + 1) * (NumDivision2 + 1) + j + 2;
                    connectivities.push_back(Node1);
                    connectivities.push_back(Node2);
                    connectivities.push_back(Node4);
                    connectivities.push_back(Node3);
                    connectivities_list.push_back(connectivities);
                }
            }
        }
        else if(mDim == 3)
        {
            int cnt = 0;
            Vector W;
            for(std::size_t i = 0; i < SamplingKnots1.size(); ++i)
            {
                for(std::size_t j = 0; j < SamplingKnots2.size(); ++j)
                {
                    for(std::size_t k = 0; k < SamplingKnots3.size(); ++k)
                    {
                        // check if this point lie in which cell
                        bool found = false;
                        for(std::vector<cell_t>::iterator it_cell = active_cells.begin(); it_cell != active_cells.end(); ++it_cell)
                        {
                            if((*it_cell)->IsCoverage(SamplingKnots1[i], SamplingKnots2[j], SamplingKnots3[k]))
                            {
                                found = true;

                                // compute the local coordinates of the sampling point
                                double local_xi = (SamplingKnots1[i] - (*it_cell)->LeftValue()) / ((*it_cell)->RightValue() - (*it_cell)->LeftValue());
                                double local_eta = (SamplingKnots2[j] - (*it_cell)->DownValue()) / ((*it_cell)->UpValue() - (*it_cell)->DownValue());
                                double local_zeta = (SamplingKnots3[k] - (*it_cell)->BelowValue()) / ((*it_cell)->AboveValue() - (*it_cell)->BelowValue());
//                                KRATOS_WATCH(local_xi)
//                                KRATOS_WATCH(local_eta)
//                                KRATOS_WATCH(local_zeta)

                                // compute the Bernstein basis function on the local coordinates
                                Vector B((mOrder1 + 1) * (mOrder2 + 1) * (mOrder3 + 1));
                                for(unsigned int i1 = 0; i1 < mOrder1 + 1; ++i1)
                                    for(unsigned int i2 = 0; i2 < mOrder2 + 1; ++i2)
                                        for(unsigned int i3 = 0; i3 < mOrder3 + 1; ++i3)
                                        {
                                            unsigned int num = (i1 * (mOrder2 + 1) + i2) * (mOrder3 + 1) + i3;
                                            double B1 = BezierUtils::bernstein2(i1, mOrder1, local_xi);
                                            double B2 = BezierUtils::bernstein2(i2, mOrder2, local_eta);
                                            double B3 = BezierUtils::bernstein2(i3, mOrder3, local_zeta);
                                            B(num) = B1 * B2 * B3;
                                        }

                                // get the extraction operator on this cell
                                Matrix C = (*it_cell)->GetExtractionOperator();

                                // compute the B-splines shape function values at the local coordinates
                                Vector N((*it_cell)->NumberOfAnchors());
                                noalias(N) = prod(C, B);

                                // compute the NURBS shape function at the local coordinates
                                (*it_cell)->GetAnchorWeights(W);
                                double Denom = inner_prod(W, N);
                                Vector R((*it_cell)->NumberOfAnchors());
                                for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                    R(k) = W(k) * N(k) / Denom;

                                // get the list of supported basis function of this cell
                                const std::vector<unsigned int>& bfs_id = (*it_cell)->GetSupportedAnchors();

                                // compute the glocal coordinates at the local coordinates
                                double X = 0.0, Y = 0.0, Z = 0.0;
                                for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                {
                                    X += R(k) * active_bfs[bfs_id[k]]->X0();
                                    Y += R(k) * active_bfs[bfs_id[k]]->Y0();
                                    Z += R(k) * active_bfs[bfs_id[k]]->Z0();
                                }

                                // add point to list
                                point_list.push_back(++cnt);
                                X_list[cnt] = X;
                                Y_list[cnt] = Y;
                                Z_list[cnt] = Z;

                                xi_list[cnt] = local_xi;
                                eta_list[cnt] = local_eta;
                                zeta_list[cnt] = 0.0;

                                parent_element_list[cnt] = (*it_cell)->Id();
                                break;
                            }
                        }
                        if(found == false)
                            KRATOS_THROW_ERROR(std::logic_error, "The sampling points is not detected in the hierarchical mesh", "")
                    }
                }
            }

            // generate the element connectivities
            std::vector<unsigned int> connectivities(8);
            for(std::size_t i = 0; i < NumDivision1; ++i)
            {
                for(std::size_t j = 0; j < NumDivision2; ++j)
                {
                    for(std::size_t k = 0; k < NumDivision3; ++k)
                    {
                        unsigned int Node1 = (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node2 = (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node3 = ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node4 = ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node5 = Node1 + 1;
                        unsigned int Node6 = Node2 + 1;
                        unsigned int Node7 = Node3 + 1;
                        unsigned int Node8 = Node4 + 1;

                        connectivities.push_back(Node1);
                        connectivities.push_back(Node2);
                        connectivities.push_back(Node4);
                        connectivities.push_back(Node3);
                        connectivities.push_back(Node5);
                        connectivities.push_back(Node6);
                        connectivities.push_back(Node8);
                        connectivities.push_back(Node7);

                        connectivities_list.push_back(connectivities);
                    }
                }
            }
        }
        
        // export to post MDPA
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for post processing of Bezier-based discretization\n";
        outfile << "//(c) 2015 Hoang Giang Bui, Ruhr-University Bochum\n";

        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        outfile << "//This file is created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900)
                << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl << std::endl;

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " " << X_list[point_list[i]]
                                     << " " << Y_list[point_list[i]]
                                     << " " << Z_list[point_list[i]]
                                     << "\n";
        outfile << "End Nodes\n\n";

        // write elements
        if(mDim == 2)
            outfile << "Begin Elements KinematicLinear2D4N\n";
        else if(mDim == 3)
            outfile << "Begin Elements KinematicLinear3D8N\n";
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")
        for(std::size_t i = 0; i < connectivities_list.size(); ++i)
        {
            outfile << (i + 1) << " 1";
            for(std::size_t j = 0; j < connectivities_list[i].size(); ++j)
                outfile << " " << connectivities_list[i][j];
            outfile << std::endl;
        }
        outfile << "End Elements\n\n";
        
        // write nodal data
        outfile << "Begin NodalData LOCAL_COORDINATES\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 [3] (" << xi_list[point_list[i]]
                                     << "," << eta_list[point_list[i]]
                                     << "," << zeta_list[point_list[i]] << ")\n";
        outfile << "End NodalData\n\n";

        outfile << "Begin NodalData PARENT_ELEMENT_ID\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 " << parent_element_list[point_list[i]] << "\n";
        outfile << "End NodalData\n\n";
        
        std::cout << "Export post MDPA to " << fn << " completed" << std::endl;
    }

    void HnMesh::FindBfs(unsigned int level, double xi_min, double xi_max, double eta_min, double eta_max, double zeta_min, double zeta_max) const
    {
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
        {
            if((*it_lvl)->Id() == level)
            {
                std::cout << "Active bfs in level " << level << " in ([" << xi_min << " " << xi_max << "],[" << eta_min << " " << eta_max << "],[" << zeta_min << " " << zeta_max << "]):";
                for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
                {
                    if((*it_bf)->IsActive())
                    {
                        double bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max;
                        (*it_bf)->GetBoundingBox(bb_xi_min, bb_xi_max, bb_eta_min, bb_eta_max, bb_zeta_min, bb_zeta_max);
                        if( bb_xi_min >= xi_min && bb_xi_max <= xi_max
                         && bb_eta_min >= eta_min && bb_eta_max <= eta_max
                         && bb_zeta_min >= zeta_min && bb_zeta_max <= zeta_max )
                            std::cout << " " << (*it_bf)->Id();
                    }
                }
                std::cout << std::endl;
                break;
            }
        }
    }

    void HnMesh::PrintBfChildren(unsigned int Id) const
    {
        std::cout << "Children of bf " << Id << ":";
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
                if((*it_bf)->Id() == Id)
                {
                    for(HnBasisFunction::bf_const_iterator it_bf_child = (*it_bf)->bf_const_begin(); it_bf_child != (*it_bf)->bf_const_end(); ++it_bf_child)
                        std::cout << " " << (*it_bf_child)->Id();
                    break;
                }
        std::cout << std::endl;
    }

    void HnMesh::PrintLocalKnotVector(unsigned int Id) const
    {
        std::cout << "Local knot vector of bf " << Id << ":" << std::endl;
        for(level_const_iterator it_lvl = level_const_begin(); it_lvl != level_const_end(); ++it_lvl)
            for(HnLevel::bf_const_iterator it_bf = (*it_lvl)->bf_const_begin(); it_bf != (*it_lvl)->bf_const_end(); ++it_bf)
                if((*it_bf)->Id() == Id)
                {
                    Vector knots1;
                    (*it_bf)->GetLocalKnots(1, knots1);
                    KRATOS_WATCH(knots1)

                    Vector knots2;
                    (*it_bf)->GetLocalKnots(2, knots2);
                    KRATOS_WATCH(knots2)

                    Vector knots3;
                    (*it_bf)->GetLocalKnots(3, knots3);
                    KRATOS_WATCH(knots3)

                    break;
                }
        std::cout << std::endl;
    }

} // end namespace Kratos

#undef MDPA_NODE_RENUMBERING
#undef MDPA_CELL_RENUMBERING
#undef ENABLE_PROFILING
#undef DEBUG_REFINE


#include "hn_basis_function.h"

// #define DEBUG_BEZIER_EXTRACTION

namespace Kratos
{
    void HnBasisFunction::ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2)
    {
        // a priori check
        for(std::size_t i = 0; i < mLocalKnots1.size(); ++i)
            if(mLocalKnots1[i] > p_cell->LeftValue() && mLocalKnots1[i] < p_cell->RightValue())
            {
                std::stringstream ss;
                ss << "Error: the cell is not contained in one knot span in u-direction of the basis function" << std::endl;
                ss << "LeftValue: " << p_cell->LeftValue() << ", RightValue: " << p_cell->RightValue() << std::endl;
                ss << "mLocalKnots1:";
                for(std::size_t j = 0; j < mLocalKnots1.size(); ++j)
                    ss << " " << mLocalKnots1[j];
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        for(std::size_t i = 0; i < mLocalKnots2.size(); ++i)
            if(mLocalKnots2[i] > p_cell->DownValue() && mLocalKnots2[i] < p_cell->UpValue())
            {
                std::stringstream ss;
                ss << "Error: the cell is not contained in one knot span in v-direction of the basis function" << std::endl;
                ss << "DownValue: " << p_cell->DownValue() << ", UpValue: " << p_cell->UpValue() << std::endl;
                ss << "mLocalKnots2:";
                for(std::size_t j = 0; j < mLocalKnots2.size(); ++j)
                    ss << " " << mLocalKnots2[j];
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

        #ifdef DEBUG_BEZIER_EXTRACTION
        std::cout << "Bezier extraction debug at bf " << Id() << ":" << std::endl;
        std::cout << "mLocalKnots1:";
        for(std::size_t i = 0; i < mLocalKnots1.size(); ++i)
            std::cout << " " << mLocalKnots1[i];
        std::cout << std::endl;
        std::cout << "mLocalKnots2:";
        for(std::size_t i = 0; i < mLocalKnots2.size(); ++i)
            std::cout << " " << mLocalKnots2[i];
        std::cout << std::endl;
        KRATOS_WATCH(p_cell->LeftValue())
        KRATOS_WATCH(p_cell->RightValue())
        KRATOS_WATCH(p_cell->DownValue())
        KRATOS_WATCH(p_cell->UpValue())
        #endif

        // compute the inserted knot vector
        std::vector<double> ins_knots1;
//        std::vector<int> ins_span1;
        for(std::size_t i = 0; i < mLocalKnots1.size() - 1; ++i)
            if(p_cell->LeftValue() > mLocalKnots1[i] && p_cell->LeftValue() < mLocalKnots1[i+1])
            {
                ins_knots1.push_back(p_cell->LeftValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
                break;
            }
        for(std::size_t i = 0; i < mLocalKnots1.size() - 1; ++i)
            if(p_cell->RightValue() > mLocalKnots1[i] && p_cell->RightValue() < mLocalKnots1[i+1])
            {
                ins_knots1.push_back(p_cell->RightValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
                break;
            }

        std::vector<double> ins_knots2;
//        std::vector<int> ins_span2;
        for(std::size_t i = 0; i < mLocalKnots2.size() - 1; ++i)
            if(p_cell->DownValue() > mLocalKnots2[i] && p_cell->DownValue() < mLocalKnots2[i+1])
            {
                ins_knots2.push_back(p_cell->DownValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
                break;
            }
        for(std::size_t i = 0; i < mLocalKnots2.size() - 1; ++i)
            if(p_cell->UpValue() > mLocalKnots2[i] && p_cell->UpValue() < mLocalKnots2[i+1])
            {
                ins_knots2.push_back(p_cell->UpValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
                break;
            }

        #ifdef DEBUG_BEZIER_EXTRACTION
        std::cout << "ins_knots1:";
        for(std::size_t i = 0; i < ins_knots1.size(); ++i)
            std::cout << " " << ins_knots1[i];
        std::cout << std::endl;
//        std::cout << "ins_span1:";
//        for(std::size_t i = 0; i < ins_span1.size(); ++i)
//            std::cout << " " << ins_span1[i];
//        std::cout << std::endl;

        std::cout << "ins_knots2:";
        for(std::size_t i = 0; i < ins_knots2.size(); ++i)
            std::cout << " " << ins_knots2[i];
        std::cout << std::endl;
//        std::cout << "ins_span2:";
//        for(std::size_t i = 0; i < ins_span2.size(); ++i)
//            std::cout << " " << ins_span2[i];
//        std::cout << std::endl;
        #endif

        // compute the Bezier extraction operator        
        std::vector<Vector> Crows;
        int nb_xi, nb_eta;
        Vector Ubar_xi, Ubar_eta;
//        BezierUtils::bezier_extraction_tsplines_2d(Crows,
//                                                   nb_xi,
//                                                   nb_eta,
//                                                   Ubar_xi,
//                                                   Ubar_eta,
//                                                   mLocalKnots1,
//                                                   mLocalKnots2,
//                                                   ins_knots1,
//                                                   ins_knots2,
//                                                   ins_span1,
//                                                   ins_span2,
//                                                   p1,
//                                                   p2);

        BezierUtils::bezier_extraction_local_2d(Crows,
                                                nb_xi,
                                                nb_eta,
                                                Ubar_xi,
                                                Ubar_eta,
                                                mLocalKnots1,
                                                mLocalKnots2,
                                                ins_knots1,
                                                ins_knots2,
                                                p1,
                                                p2);

        #ifdef DEBUG_BEZIER_EXTRACTION
        std::cout << "Crows:" << std::endl;
        for(std::size_t i = 0; i < Crows.size(); ++i)
            std::cout << Crows[i] << std::endl;
        #endif

        // extract the correct row
        int span1;// knot span in u-direction of the cell w.r.t basis function support
        int span2; // knot span in v-direction of the cell w.r.t basis function support
        #ifdef DEBUG_BEZIER_EXTRACTION
        KRATOS_WATCH(Ubar_xi)
        KRATOS_WATCH(Ubar_eta)
        #endif
        std::set<double> Ubar_xi_unique(Ubar_xi.begin(), Ubar_xi.end());
        std::set<double> Ubar_eta_unique(Ubar_eta.begin(), Ubar_eta.end());
        std::vector<double> Ubar_xi_unique_vector(Ubar_xi_unique.begin(), Ubar_xi_unique.end());
        std::vector<double> Ubar_eta_unique_vector(Ubar_eta_unique.begin(), Ubar_eta_unique.end());
        span1 = this->FindSpanLocal(p_cell->LeftValue(), Ubar_xi_unique_vector) - 1;
        span2 = this->FindSpanLocal(p_cell->DownValue(), Ubar_eta_unique_vector) - 1;

        #ifdef DEBUG_BEZIER_EXTRACTION
        KRATOS_WATCH(span1)
        KRATOS_WATCH(span2)
        #endif
//        int span = span2 * nb_xi + span1;
        int span = span1 * nb_eta + span2;
        Crow.resize(Crows[span].size());
        #ifdef DEBUG_BEZIER_EXTRACTION
        KRATOS_WATCH(span)
        KRATOS_WATCH(Crows.size())
        #endif
        std::copy(Crows[span].begin(), Crows[span].end(), Crow.begin());
        #ifdef DEBUG_BEZIER_EXTRACTION
        std::cout << "----------------------" << std::endl;
        #endif
    }
    
    void HnBasisFunction::ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2, int p3)
    {
        // TODO
    }
}

#undef DEBUG_BEZIER_EXTRACTION


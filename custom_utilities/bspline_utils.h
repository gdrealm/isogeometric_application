//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-08-18 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BSPLINE_UTILS_H_INCLUDED )
#define  KRATOS_BSPLINE_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream> 

// External includes 
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/isogeometric_math_utils.h"

#define FUNCTIONALITY_CHECK // this macro is used to enable the compulsory functionality check of the program. If we are sure if it works, then we can disable it.

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class BSplineUtils
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;
    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;

    /// Pointer definition of BSplineUtils
    KRATOS_CLASS_POINTER_DEFINITION(BSplineUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BSplineUtils()
    {}

    /// Destructor.
    virtual ~BSplineUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // original implementation by hbui
//    template<class ValuesContainerType>
//      static int FindSpan(
//        const int rN,
//        const int rP,
//        const double rXi,
//        const ValueContainerType& rU
//        )
//      {
//        int size = rU.size();
//
//        if(size != rN + rP + 1)
//        {
//            KRATOS_THROW_ERROR(std::logic_error, "Size of rU must be n + p + 1", __FUNCTION__);
//        }
//
//        double tol = 1e-12;
//
//        if(fabs(rXi - 1.0) < tol)
//            return rN - 1;
//
//        int span_a = 0;
//        int span_b = rN + rP;
//        int span_mid;
//
//        while( abs(span_a - span_b) != 1 )
//        {
//            span_mid = (span_a + span_b) / 2;
//
//            if(rXi < rU[span_mid])
//            {
//                span_b = span_mid;
//            }
//            else
//            {
//                span_a = span_mid;
//            }
//        }
//
//        return span_a;
//      }

    //implementation from the NURBS book
    // Note:
    // Algorithm A2.1 from 'The NURBS BOOK' pg68
    // as that algorithm only works for nonperiodic
    // knot vectors, nonetheless the results should 
    // be EXACTLY the same if U is nonperiodic
    // sample: rN = 4, rP = 2, rU = [0 0 0 0.5 1 1 1], rXi = 0.1 (rU.size() == rN + rP + 1)
    // note: this only works with open knot vector (i.e. 0  and 1 is repeated p + 1 times)
    template<class ValuesContainerType>
    static int FindSpan(
            const int rN,
            const int rP,
            const double rXi,
            const ValuesContainerType& rU
    )
    {
        if(rXi == rU[rN]) return rN - 1;

        int low = rP;
        int high = rN;
        int mid = (low + high) / 2;

        while( rXi < rU[mid] || rXi >= rU[mid+1] )
        {
            if(rXi < rU[mid])
            {
                high = mid;
            }
            else
            {
                low = mid;
            }
            mid = (low + high) / 2;
        }

        return mid;
    }
    
    // implementation in GeoPde, low_level_functions.cc
    // Note: this implementation has linear, rather than log complexity
//    template<class ValuesContainerType>
//    static int FindSpan(
//            const int rN,
//            const int rP,
//            const double rXi,
//            const ValuesContainerType& rU
//    )
//    {
//        int ret = 0;
//        KRATOS_WATCH(rN)
//        while ((ret++ < rN) && (rU[ret] <= rXi))
//        {
//            KRATOS_WATCH(ret)
//        }
//        if(ret == rN+1) --ret;
//        std::cout << "--------" << std::endl;
//        return (ret - 1);
//    }
  
    // Remark: this function only works correctly when rXi is in the knot span rI.
//    static void BasisFuns(
//            ValueContainerType& rS,
//            const int rI,
//            const double rXi,
//            const int rP,
//            const ValueContainerType& rU
//    )
//    {
//        rS[0] = 1.0;
//
//        for(unsigned int j = 1; j < rP + 1; ++j)
//        {
//            rS[j] = (rXi - rU[rI]) / (rU[rI + j] - rU[rI]) * rS[j-1];
//
//            for(unsigned int i = j - 1; i > 0; --i)
//            {
//                int t = j - i;
//                rS[i] = (rXi - rU[rI - t]) / (rU[rI - t + j] - rU[rI - t]) * rS[i - 1]
//                + (rU[rI - t + j + 1] - rXi) / (rU[rI - t + j + 1] - rU[rI - t + 1]) * rS[i];
//            }
//
//            rS[0] *= (rU[rI + 1] - rXi) / (rU[rI + 1] - rU[rI - j + 1]);
//        }
//    }
    //time:
    //N = 100: 0.000123978
    //N = 1000: 0.00120997
    //N = 10000: 0.013212
    //N = 100000: 0.119535
    //N = 1000000: 0.96382
    //N = 10000000: 8.45895
    //N = 100000000: 81.9643, 81.9494

    // This function works slightly faster than the above implementation
    template<class ValuesContainerType>
    static void BasisFuns(ValuesContainerType& rS,
                          const int rI,
                          const double rXi,
                          const int rP,
                          const ValuesContainerType& rU)
    {
        unsigned int j, r;
        double saved, temp;

        double left[rP + 1];
        double right[rP + 1];
        
        std::fill(left, left + rP + 1, 0.0);
        std::fill(right, right + rP + 1, 0.0);

        rS[0] = 1.0;
        for (j = 1; j <= rP; ++j)
        {
            left[j] = rXi - rU[rI + 1 - j];
            right[j] = rU[rI + j] - rXi;
            saved = 0.0;

            for (r = 0; r < j; ++r)
            {
                temp = rS[r] / (right[r + 1] + left[j - r]);
                rS[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }

            rS[j] = saved;
        }
    }
    //N = 10000000: 8.18893
    //N = 100000000: 76.7523, 81.9167

    /**
     * Computes b-spline function derivatives
     */
    template<class ValuesArrayContainerType, class ValuesContainerType>
    static void BasisFunsDer(ValuesArrayContainerType& rS,
                             int rI,
                             const double rXi,
                             const int rP,
                             const ValuesContainerType& rU,
                             const int rD)
    {
        int j, r, k, rk, pk, j1, j2, s1, s2;
        double d, temp, saved;

//        ders = zeros(nders+1,pl+1);
        noalias(rS) = ZeroMatrix(rD + 1, rP + 1);

//        ndu = zeros(pl+1,pl+1);
        Matrix ndu = ZeroMatrix(rP + 1, rP + 1);
//        left = zeros(pl+1);
        Vector left = ZeroVector(rP + 1);
//        right = zeros(pl+1);
        Vector right = ZeroVector(rP + 1);

//        a = zeros(2,pl+1);
        Matrix a = ZeroMatrix(2, rP + 1);

//        ndu(1,1) = 1;
        ndu(0, 0) = 1.0;
        rI += 1;

//        for j = 1:pl
        for (j = 1; j <= rP; ++j)
        {
//            left(j+1) = u - u_knotl(i+1-j);
            left(j) = rXi - rU(rI - j);
//            right(j+1) = u_knotl(i+j) - u;
            right(j) = rU(rI + j - 1) - rXi;
//            saved = 0;
            saved = 0.0;
//            for r = 0:j-1
            for (r = 0; r <= j - 1; ++r)
            {
//                ndu(j+1,r+1) = right(r+2) + left(j-r+1);
                ndu(j, r) = right(r + 1) + left(j - r);
//                temp = ndu(r+1,j)/ndu(j+1,r+1);
                temp = ndu(r, j - 1) / ndu(j, r);
//                ndu(r+1,j+1) = saved + right(r+2)*temp;
                ndu(r, j) = saved + right(r + 1) * temp;
//                saved = left(j-r+1)*temp;
                saved = left(j - r) * temp;
//            end
            }
//            ndu(j+1,j+1) = saved;
            ndu(j, j) = saved;
//        end
        }

//        for j = 0:pl
        for (j = 0; j <= rP; ++j)
        {
//            ders(1,j+1) = ndu(j+1,pl+1);
            rS(0, j) = ndu(j, rP);
//        end
        }

//        for r = 0:pl
        for (r = 0; r <= rP; ++r)
        {
            s1 = 0;
            s2 = 1;
//            a(1,1) = 1;
            a(0, 0) = 1.0;
//            for k = 1:nders //compute kth derivative
            for (k = 1; k <= rD; ++k)
            {
                d = 0.0;
                rk = r - k;
                pk = rP - k;
                if(r >= k)
                {
//                    a(s2+1 , 1) = a(s1+1,1)/ndu(pk+2,rk+1);
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
//                    d = a(s2+1,1)*ndu(rk+1,pk+1);
                    d = a(s2, 0) * ndu(rk, pk);
                }
                if(rk >= -1)
                {
                    j1 = 1;
                }
                else
                {
                    j1 = -rk;
                }
                if((r-1) <= pk)
                {
                    j2 = k - 1;
                }
                else
                {
                    j2 = rP - r;
                }
//                for j = j1:j2
                for (j = j1; j <= j2; ++j)
                {
//                    a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
//                    d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
                    d = d + a(s2, j) * ndu(rk + j, pk);
//                end
                }
                if (r <= pk)
                {
//                    a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                    a(s2, k) = -a(s1, k - 1) / ndu(pk+1, r);
//                    d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
                    d = d + a(s2, k) * ndu(r, pk);
//                end
                }
//                ders(k+1,r+1) = d;
                rS(k, r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
//            end
            }
//        end
        }

        r = rP;
//        for k = 1:nders
        for (k = 1; k <= rD; ++k)
        {
//            for j = 0:pl
            for (j = 0; j <= rP; ++j)
            {
//                ders(k+1,j+1) = ders(k+1,j+1)*r;
                rS(k, j) = rS(k, j) * r;
//            end
            }
            r = r * (rP - k);
//        end
        }

    }

    /// Compute the refinement coefficients from knot insertion refinement for 1D
    template<class MatrixType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients1D(MatrixType& D,
                                                           ValuesContainerType& new_knots,
                                                           int p,
                                                           const ValuesContainerType& knots,
                                                           double k)
    {
        // compute the number of basis function
        int n = knots.size() - p - 1;

        // find the span of the inserted knot
        int s = FindSpan(n, p, k, knots);
//        std::cout << "knots:";
//        for(int i = 0; i < knots.size(); ++i)
//            std::cout << " " << knots[i];
//        std::cout << std::endl;
//        std::cout << "k=" << k << ",s=" << s << std::endl;
        
        // form the new knot vector
        new_knots.resize(knots.size() + 1);
        for(int i = 0; i < s+1; ++i)
            new_knots[i] = knots[i];
        new_knots[s+1] = k;
        for(int i = s+2; i < knots.size() + 1; ++i)
            new_knots[i] = knots[i-1];

//        std::cout << "new_knots:";
//        for(int i = 0; i < new_knots.size(); ++i)
//            std::cout << " " << new_knots[i];
//        std::cout << std::endl;

        // initialize and compute coefficient matrix
        D.resize(n, n+1);
        noalias(D) = ZeroMatrix(n, n+1);
        for(int i = 0; i < s-p; ++i)
            D(i, i) = 1.0;
        for(int i = s-p; i < s+1; ++i)
        {
            D(i, i) = (k - new_knots[i]) / (new_knots[i+p+1] - new_knots[i]);
            D(i, i+1) = (new_knots[i+p+2] - k) / (new_knots[i+p+2] - new_knots[i+1]);
        }
        for(int i = s+1; i < n; ++i)
            D(i, i+1) = 1.0;
//        KRATOS_WATCH(D)
//        std::cout << "-----------------" << std::endl;
    }

    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2>
    static void ComputeBsplinesKnotInsertionCoefficients1D(MatrixType& D,
                                                           ValuesContainerType& new_knots,
                                                           int p,
                                                           const ValuesContainerType& knots,
                                                           const ValuesContainerType2& ins_knots)
    {
        // compute the number of basis function
        int n = knots.size() - p - 1;

        // initialize coefficient matrix
        D.resize(n, n);
        noalias(D) = IdentityMatrix(n, n);

        // copy the knot vector
        new_knots.resize(knots.size());
        std::copy(knots.begin(), knots.end(), new_knots.begin());

        ValuesContainerType tmp_knots;
        tmp_knots.resize(new_knots.size());
        std::copy(new_knots.begin(), new_knots.end(), tmp_knots.begin());

        // insert individual knots and concatenate coefficient matrix
        MatrixType Di, Dt;
        for(std::size_t i = 0; i < ins_knots.size(); ++i)
        {
            ComputeBsplinesKnotInsertionCoefficients1D(Di, new_knots, p, tmp_knots, ins_knots[i]);

            Dt.resize(D.size1(), Di.size2());
            noalias(Dt) = prod(D, Di);

            D.resize(Dt.size1(), Dt.size2());
            noalias(D) = Dt;

            tmp_knots.resize(new_knots.size());
            std::copy(new_knots.begin(), new_knots.end(), tmp_knots.begin());
        }
//        KRATOS_WATCH(D)
    }

    template<class MatrixType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients2D(MatrixType& D,
                                                           ValuesContainerType& new_knots1,
                                                           ValuesContainerType& new_knots2,
                                                           int p1,
                                                           int p2,
                                                           const ValuesContainerType& knots1,
                                                           const ValuesContainerType& knots2,
                                                           const ValuesContainerType& ins_knots1,
                                                           const ValuesContainerType& ins_knots2)
    {
        MatrixType D1, D2;
        ComputeBsplinesKnotInsertionCoefficients1D(D1, new_knots1, p1, knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1D(D2, new_knots2, p2, knots2, ins_knots2);
//        KRATOS_WATCH(D1)
//        KRATOS_WATCH(D2)
        IsogeometricMathUtils::outer_prod_mat(D, D2, D1);
//        KRATOS_WATCH(D)
    }

    template<class MatrixType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients3D(MatrixType& D,
                                                           ValuesContainerType& new_knots1,
                                                           ValuesContainerType& new_knots2,
                                                           ValuesContainerType& new_knots3,
                                                           int p1,
                                                           int p2,
                                                           int p3,
                                                           const ValuesContainerType& knots1,
                                                           const ValuesContainerType& knots2,
                                                           const ValuesContainerType& knots3,
                                                           const ValuesContainerType& ins_knots1,
                                                           const ValuesContainerType& ins_knots2,
                                                           const ValuesContainerType& ins_knots3)
    {
        MatrixType D1, D2, D3;
        ComputeBsplinesKnotInsertionCoefficients1D(D1, new_knots1, p1, knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1D(D2, new_knots2, p2, knots2, ins_knots2);
        ComputeBsplinesKnotInsertionCoefficients1D(D3, new_knots3, p3, knots3, ins_knots3);
        
        MatrixType tmp;
        IsogeometricMathUtils::outer_prod_mat(tmp, D2, D1);
        IsogeometricMathUtils::outer_prod_mat(D, D3, tmp);
    }

    template<class MatrixType, class ValuesContainerType>
    static void ComputeNURBSKnotInsertionCoefficients1D(MatrixType& D,
                                                        ValuesContainerType& new_knots,
                                                        ValuesContainerType& new_weights,
                                                        int p,
                                                        const ValuesContainerType& knots,
                                                        const ValuesContainerType& ins_knots,
                                                        const ValuesContainerType& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients1D(D, new_knots, p, knots, ins_knots);

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
    }

    template<class MatrixType, class ValuesContainerType>
    static void ComputeNURBSKnotInsertionCoefficients2D(MatrixType& D,
                                                        ValuesContainerType& new_knots1,
                                                        ValuesContainerType& new_knots2,
                                                        ValuesContainerType& new_weights,
                                                        int p1,
                                                        int p2,
                                                        const ValuesContainerType& knots1,
                                                        const ValuesContainerType& knots2,
                                                        const ValuesContainerType& ins_knots1,
                                                        const ValuesContainerType& ins_knots2,
                                                        const ValuesContainerType& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients2D(D, new_knots1, new_knots2, p1, p2, knots1, knots2, ins_knots1, ins_knots2);
//        std::cout << "B-splines D:" << D << std::endl;

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }
//        std::cout << "new_weights:";
//        for(int i = 0; i < new_weights.size(); ++i)
//            std::cout << " " << new_weights[i];
//        std::cout << std::endl;

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
//        KRATOS_WATCH(D)
    }

    template<class MatrixType, class ValuesContainerType>
    static void ComputeNURBSKnotInsertionCoefficients3D(MatrixType& D,
                                                        ValuesContainerType& new_knots1,
                                                        ValuesContainerType& new_knots2,
                                                        ValuesContainerType& new_knots3,
                                                        ValuesContainerType& new_weights,
                                                        int p1,
                                                        int p2,
                                                        int p3,
                                                        const ValuesContainerType& knots1,
                                                        const ValuesContainerType& knots2,
                                                        const ValuesContainerType& knots3,
                                                        const ValuesContainerType& ins_knots1,
                                                        const ValuesContainerType& ins_knots2,
                                                        const ValuesContainerType& ins_knots3,
                                                        const ValuesContainerType& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients3D(D, new_knots1, new_knots2, new_knots3, p1, p2, p3, knots1, knots2, knots3, ins_knots1, ins_knots2, ins_knots3);

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
    }

    /**
        TODO: to be tested
     */
    template<class VectorType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients1DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots,
                                                                int p,
                                                                const ValuesContainerType& local_knots,
                                                                const ValuesContainerType& ins_knots)
    {
        // compute the extended knot vector
        ValuesContainerType Ubar;
        int nt;
        IsogeometricMathUtils::compute_extended_knot_vector(Ubar, nt, local_knots, p);
//        KRATOS_WATCH(nt)

        // compute the refinement matrix for extended knot vector
        Matrix M;
        ComputeBsplinesKnotInsertionCoefficients1D(M, new_knots, p, Ubar, ins_knots);
//        KRATOS_WATCH(M)

        // extract the correct row from the refinement matrix
        int num_basis = ins_knots.size() + 1;
        if(D.size() != num_basis)
            D.resize(num_basis);
        for(unsigned int i = 0; i < num_basis; ++i)
            D[i] = M(nt, nt + i);

        #ifdef FUNCTIONALITY_CHECK
        for(unsigned int i = 0; i < nt; ++i)
            if(M(nt, i) != 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "M(nt,..) is nonzero at", i)
        for(unsigned int i = 0; i < num_basis; ++i)
            if(D(i) == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "D(..) is zero at", i)
        for(unsigned int i = nt + num_basis; i < M.size2(); ++i)
            if(M(nt, i) != 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "M(nt,..) is nonzero at", i)
        #endif
    }

    template<class VectorType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients2DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots1,
                                                                ValuesContainerType& new_knots2,
                                                                int p1,
                                                                int p2,
                                                                const ValuesContainerType& local_knots1,
                                                                const ValuesContainerType& local_knots2,
                                                                const ValuesContainerType& ins_knots1,
                                                                const ValuesContainerType& ins_knots2)
    {
        VectorType D1, D2;
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D1, new_knots1, p1, local_knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D2, new_knots2, p2, local_knots2, ins_knots2);
        IsogeometricMathUtils::outer_prod_vec(D, D2, D1);
    }

    template<class VectorType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients3DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots1,
                                                                ValuesContainerType& new_knots2,
                                                                ValuesContainerType& new_knots3,
                                                                int p1,
                                                                int p2,
                                                                int p3,
                                                                const ValuesContainerType& local_knots1,
                                                                const ValuesContainerType& local_knots2,
                                                                const ValuesContainerType& local_knots3,
                                                                const ValuesContainerType& ins_knots1,
                                                                const ValuesContainerType& ins_knots2,
                                                                const ValuesContainerType& ins_knots3)
    {
        VectorType D1, D2, D3;
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D1, new_knots1, p1, local_knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D2, new_knots2, p2, local_knots2, ins_knots2);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D3, new_knots3, p3, local_knots3, ins_knots3);

        VectorType aux;
        IsogeometricMathUtils::outer_prod_vec(aux, D2, D1);
        IsogeometricMathUtils::outer_prod_vec(D, D3, aux);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Testing
    ///@{

    void test_ComputeBsplinesKnotInsertionCoefficients1DLocal();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "BSplineUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplineUtils";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BSplineUtils& operator=(BSplineUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BSplineUtils(BSplineUtils const& rOther)
    {
    }

    ///@}

}; // Class BSplineUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BSplineUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const BSplineUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef FUNCTIONALITY_CHECK

#endif // KRATOS_BSPLINE_UTILS_H_INCLUDED  defined

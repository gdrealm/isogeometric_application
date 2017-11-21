//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <vector>

// External includes

// Project includes
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"

namespace Kratos
{

template<>
std::vector<double> BSplinesFESpace<2>::GetValue(const std::vector<double>& xi) const
{
    // locate the knot span
    int Span[2];
    Span[0] = BSplineUtils::FindSpan(this->Number(0), this->Order(0), xi[0], this->KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(this->Number(1), this->Order(1), xi[1], this->KnotVector(1));

    // compute the non-zero shape function
    std::vector<double> ShapeFunctionValues1(this->Order(0) + 1);
    std::vector<double> ShapeFunctionValues2(this->Order(1) + 1);

    BSplineUtils::BasisFuns(ShapeFunctionValues1, Span[0], xi[0], this->Order(0), this->KnotVector(0));
    BSplineUtils::BasisFuns(ShapeFunctionValues2, Span[1], xi[1], this->Order(1), this->KnotVector(1));

    // rearrange the shape functions
    std::vector<double> results(this->TotalNumber());
    std::fill(results.begin(), results.end(), 0.0);

    int Start[2];
    Start[0] = Span[0] - this->Order(0);
    Start[1] = Span[1] - this->Order(1);

    double N1;
    double N2;

    unsigned int i, j, Index;
    for(i = Start[0]; i <= Span[0]; ++i)
    {
        for(j = Start[1]; j <= Span[1]; ++j)
        {
            Index = BSplinesIndexingUtility::Index2D(i+1, j+1, this->Number(0), this->Number(1));

            N1 = ShapeFunctionValues1[i - Start[0]];
            N2 = ShapeFunctionValues2[j - Start[1]];

            results[Index] = N1 * N2;
        }
    }

    return results;
}

template<>
std::vector<double> BSplinesFESpace<3>::GetValue(const std::vector<double>& xi) const
{
    // locate the knot span
    int Span[3];
    Span[0] = BSplineUtils::FindSpan(this->Number(0), this->Order(0), xi[0], this->KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(this->Number(1), this->Order(1), xi[1], this->KnotVector(1));
    Span[2] = BSplineUtils::FindSpan(this->Number(2), this->Order(2), xi[2], this->KnotVector(2));

    // compute the non-zero shape function
    std::vector<double> ShapeFunctionValues1(this->Order(0) + 1);
    std::vector<double> ShapeFunctionValues2(this->Order(1) + 1);
    std::vector<double> ShapeFunctionValues3(this->Order(2) + 1);

    BSplineUtils::BasisFuns(ShapeFunctionValues1, Span[0], xi[0], this->Order(0), this->KnotVector(0));
    BSplineUtils::BasisFuns(ShapeFunctionValues2, Span[1], xi[1], this->Order(1), this->KnotVector(1));
    BSplineUtils::BasisFuns(ShapeFunctionValues3, Span[2], xi[2], this->Order(2), this->KnotVector(2));

    // rearrange the shape functions
    std::vector<double> results(this->TotalNumber());
    std::fill(results.begin(), results.end(), 0.0);

    int Start[3];
    Start[0] = Span[0] - this->Order(0);
    Start[1] = Span[1] - this->Order(1);
    Start[2] = Span[2] - this->Order(2);

    double N1;
    double N2;
    double N3;

    unsigned int i, j, k, Index;
    for(i = Start[0]; i <= Span[0]; ++i)
    {
        for(j = Start[1]; j <= Span[1]; ++j)
        {
            for(k = Start[2]; k <= Span[2]; ++k)
            {
                Index = BSplinesIndexingUtility::Index3D(i+1, j+1, k+1, this->Number(0), this->Number(1), this->Number(2));

                N1 = ShapeFunctionValues1[i - Start[0]];
                N2 = ShapeFunctionValues2[j - Start[1]];
                N3 = ShapeFunctionValues3[k - Start[2]];

                results[Index] = N1 * N2 * N3;
            }
        }
    }


    return results;
}

} // namespace Kratos.


#include "nurbscurve.h"

using namespace std;
using namespace geomodcore;


NURBSCurve::NURBSCurve(const Polygon& poly, const Knots& knotValues, const unsigned int K) 
	: BSplineCurve(poly, knotValues, K) 
{
}

NURBSCurve::NURBSCurve(const PolygonBasedCurve& curve, const Knots& knotValues, const unsigned int K) 
	: BSplineCurve(curve, knotValues, K)
{
}


NURBSCurve::NURBSCurve(const BSplineCurve& curve) 
	: BSplineCurve(curve)
{
}


NURBSCurve::NURBSCurve(const Polygon& poly, const unsigned int K)
	: BSplineCurve(poly, K) 
{
	const Knots& knots = initUniformKnotValues((unsigned int) poly.size(), K);
	this->knots.assign(knots.cbegin(), knots.cend());
}

NURBSCurve::NURBSCurve(const PolygonBasedCurve& curve, const unsigned int K)
	: BSplineCurve(curve,K)
{
	const Knots& knots = initUniformKnotValues(curve.getNrOfControlPoints(), K);
	this->knots.assign(knots.cbegin(), knots.cend());
}

NURBSCurve::NURBSCurve(void) 
{
}

NURBSCurve::NURBSCurve(const NURBSCurve& curve) 
	: BSplineCurve(curve) 
{
}

NURBSCurve::~NURBSCurve(void) 
{
}

NURBSCurve& NURBSCurve::operator= (const NURBSCurve& curve) 
{
	NURBSCurve* newCurve = new NURBSCurve(curve);
	return *newCurve;
}

GPoint NURBSCurve::getPoint(double t) const 
{
	const Polygon& P = polygon;
	const size_t N = P.size();
	double p0=0, p1=0, p2=0;
	for (unsigned int i = 0; i < N; i++) {
		const double b = BSplineCurve::BSplineFunction(t, i, K, (unsigned int) N, K, knots);
		p0 += P[i][0] * b;
		p1 += P[i][1] * b;
		p2 += P[i][2] * b;
	}
	GPoint p = {p0,p1,p2};
	return p;
}

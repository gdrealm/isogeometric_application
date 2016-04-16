#include "beziercurve.h"

using namespace std;
using namespace geomodcore;
using namespace boost;

BezierCurve::BezierCurve(const Polygon& poly) 
	: PolygonBasedCurve(poly)
{
}

BezierCurve::BezierCurve(const PolygonBasedCurve& curve) 
	: PolygonBasedCurve(curve)
{
}


BezierCurve::BezierCurve(void) 
{
}

BezierCurve::BezierCurve(const BezierCurve& curve)
	: PolygonBasedCurve(curve)
{
}

BezierCurve::~BezierCurve(void) 
{
}

BezierCurve& BezierCurve::operator= (const BezierCurve& curve) 
{
	BezierCurve* newCurve = new BezierCurve(curve);
	return *newCurve;
}

GPoint BezierCurve::getPoint(double t) const 
{
	const Polygon& P = polygon;
	GPoint p = {0,0,0};
	for (unsigned int i = 0; i < polygon.size(); i++) {
		const double b = bernstein(i, polygon.size(), t);
		p[0] += P[i][0] * b;
		p[1] += P[i][1] * b;
		p[2] += P[i][2] * b;
	}
	return p;
}

void BezierCurve::increaseDegree(void)
{
	Polygon& P = polygon;
	const unsigned int n = P.size();
	P.push_back(P[n-1]);
	for (int i = n; i > 0; i--) {
		const double c = boost::math::binomial_coefficient<double>(n+1, i);
		const double b0 = boost::math::binomial_coefficient<double>(n, i-1) / c;
		const double b1 = boost::math::binomial_coefficient<double>(n, i) / c;
		GPoint p = {b0 * P[i+1][0]  + b1 * P[i][0], b0 * P[i+1][1]  + b1 * P[i][1], b0 * P[i+1][2]  + b1 * P[i][2]};
		P[i] = p;
	}
}

MinMaxP1d BezierCurve::minmaxT() const 
{
	return MinMaxP1d(0.0, 1.0);
}


BezierCurve* BezierCurve::getDerivation(void) const
{
	const Polygon& P = polygon;

	Polygon derivPoly;
	for (unsigned int i = 1; i < P.size(); i++) {
		GPoint p = { P[i][0] - P[i-1][0], P[i][1] - P[i-1][1], P[i][2] - P[i-1][2] };
		derivPoly.push_back(p);
	}
	return new BezierCurve(derivPoly);
};

double BezierCurve::getArcLength(void) const 
{
	cout << "Unsupported Operation" << endl;
	return NAN;
}
#include "lagrangecurve.h"

using namespace std;
using namespace geomodcore;


LagrangeCurve::LagrangeCurve(const Polygon& poly) 
	: PolygonBasedCurve(poly)
{
}

LagrangeCurve::LagrangeCurve(const PolygonBasedCurve& curve) 
	: PolygonBasedCurve(curve)
{
}


LagrangeCurve::LagrangeCurve(void) 
{
}

LagrangeCurve::LagrangeCurve(const LagrangeCurve& curve) 
	: PolygonBasedCurve(curve)
{
}

LagrangeCurve::~LagrangeCurve(void)
{
}

LagrangeCurve& LagrangeCurve::operator= (const LagrangeCurve& curve) 
{
	LagrangeCurve* newCurve = new LagrangeCurve(curve);
	return *newCurve;
}

GPoint LagrangeCurve::getPoint(double t) const
{
	const Polygon& P = polygon;
	double p0=0, p1=0, p2=0;
	for (unsigned int i = 0; i < P.size(); i++)
	{
		const double lVal = LagrangeCurve::lagrange(i, t, P.size());
		p0 += P[i][0] * lVal;
		p1 += P[i][1] * lVal;
		p2 += P[i][2] * lVal;
	}
	GPoint p = {p0,p1,p2};
	return p;
}


MinMaxP1d LagrangeCurve::minmaxT() const 
{
	return MinMaxP1d(0.0, 1.0);
}

double LagrangeCurve::lagrange(const unsigned int i, const double t, const unsigned int K) 
{
	double n = 1.0;
	double z = 1.0;

	const double ti = i / (double) K;

	for (unsigned int j = 0; j <= K; j++) {
		if (i == j) {
			continue;
		}

		const double tj = j / (double) K;

		n *= (t - tj);
		z *= (ti - tj);
	}
	return n / z;
}

/**
* Bestimmung der abgeleiteten Kurve.
*/
LagrangeCurve* LagrangeCurve::getDerivation(void) const
{
	cout << "Unsupported Operation" << endl;
	return nullptr;
}

/**
* Berechnet die Bogenlänge der Kurve aus der Ableitung.
*/
double LagrangeCurve::getArcLength(void) const
{
	cout << "Unsupported Operation" << endl;
	return NAN;
}
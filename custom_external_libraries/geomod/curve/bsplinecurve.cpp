#include "bsplinecurve.h"

using namespace std;
using namespace geomodcore;

BSplineCurve::BSplineCurve(const Polygon& poly, const unsigned int K) 
	: PolygonBasedCurve(poly), K(K)
{
	const Knots knots = initUniformKnotValues((unsigned int) poly.size(), K);
	this->knots.assign(knots.cbegin(), knots.cend());
}
BSplineCurve::BSplineCurve(const PolygonBasedCurve& curve, const unsigned int K) 
	: PolygonBasedCurve(curve), K(K)
{
	const Knots knots = initUniformKnotValues(curve.getNrOfControlPoints(), K);
	this->knots.assign(knots.cbegin(), knots.cend());
}


BSplineCurve::BSplineCurve(const Polygon& poly, const Knots& knots, const unsigned int K) 
	: PolygonBasedCurve(poly), K(K)
{
	this->knots.assign(knots.cbegin(), knots.cend());
	if (knots.empty()) {
		initUniformKnotValues((unsigned int) poly.size(), K);
	}
}
BSplineCurve::BSplineCurve(const PolygonBasedCurve& curve, const Knots& knots, const unsigned int K) 
	: PolygonBasedCurve(curve), K(K)
{
	this->knots.assign(knots.cbegin(), knots.cend());
	if (knots.empty()) {
		initUniformKnotValues(curve.getNrOfControlPoints(), K);
	}
}

BSplineCurve::BSplineCurve(void) {
	K = 3;
	initUniformKnotValues(0,K);
}

BSplineCurve::BSplineCurve(const BSplineCurve& curve) :
PolygonBasedCurve(curve), K(curve.getK())
{
	const Knots& knots = curve.getKnotValuesConstRef();
	this->knots.assign(knots.cbegin(), knots.cend());
}

BSplineCurve::~BSplineCurve(void) {
}

BSplineCurve& BSplineCurve::operator= (const BSplineCurve& curve) {
	BSplineCurve* newCurve = new BSplineCurve(curve);
	return *newCurve;
}

GPoint BSplineCurve::getPoint(double t) const {
	const Polygon& P = polygon;
	const size_t N = polygon.size();
	double p0 = 0,p1 = 0,p2 = 0;
	for (unsigned int i = 0; i < N; i++) {
		const double b = BSplineCurve::BSplineFunction(t, i, K, (unsigned int) N, K, knots);
		p0 += P[i][0] * b;
		p1 += P[i][1] * b;
		p2 += P[i][2] * b;
	}
	GPoint p = {p0,p1,p2};
	return p;
}

MinMaxP1d BSplineCurve::minmaxT() const {
	return MinMaxP1d(knots.front(), knots.back() - geomodcore::EPSILON);
}

unsigned int BSplineCurve::getK() const
{
	return K;
}

void BSplineCurve::setK(const unsigned int K) {
	assert((K > 0) && (K <= polygon.size()));
	this->K = K;
}

Knots BSplineCurve::initUniformKnotValues(const unsigned int n, const unsigned int K)
{
	Knots knots;
	knots.assign(K-1,0);
	for (unsigned int k = 0; k <= n + 2 - K; k++) {
		knots.push_back((double) k);
	}
	knots.insert(knots.end(), K-1, (double) n + 2 - K);
	return knots;
};


/**
* Bestimmung der abgeleiteten Kurve.
*/
BSplineCurve* BSplineCurve::getDerivation(void) const
{
	cout << "Unsupported Operation" << endl;
	return nullptr;
}

/**
* Berechnet die Bogenlänge der Kurve aus der Ableitung.
*/
double BSplineCurve::getArcLength(void) const
{
	cout << "Unsupported Operation" << endl;
	return NAN;
}

double BSplineCurve::BSplineFunction(const double u, const unsigned int i, const unsigned int r, const unsigned int N, const unsigned int K, const Knots& knots) {

	// fehlerhafter Knotenvektor
	assert(knots.size() == N + K + 1);

	// Abbruchbedingung der Rekursion, wenn stückweise konstante Funktionen erreicht werden
	if (r == 0) {
		if (u >= knots[i] && u < knots[i+1]) {
			return 1.0;
		} else {
			return 0.0;
		}
	} else {
		// Rekursionsschritt
		const double tmp0 = knots[i+r] - knots[i];
		const double tmp1 = knots[i+r+1] - knots[i+1];

		const double alpha1 = (tmp0 != 0.0 ? (u - knots[i]) / tmp0 : 0.0);
		const double alpha2 = (tmp1 != 0.0 ? (knots[i+r+1] - u) / tmp1 : 0.0);

		const double B1 = (alpha1 == 0.0 ? 0.0 : BSplineFunction(u, i, r - 1, N, K, knots));
		const double B2 = (alpha2 == 0.0 ? 0.0 : BSplineFunction(u, i + 1, r - 1, N, K, knots));
		return alpha1 * B1 + alpha2 * B2;
	}
}

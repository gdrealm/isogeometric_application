#include "parametriccurve.h"

using namespace std;
using namespace geomodcore;


Polygon* ParametricCurve::getDiscretization(const unsigned int n) const {
	assert (n > 0);
	MinMaxP1d minmax = minmaxT();
	const double dt = (minmax.second - minmax.first) / n;
	Polygon* poly = new Polygon();
	for (unsigned int i = 0; i <= n; i++) {
		poly->push_back(getPoint(i * dt));
	}
	return poly;
}

double ParametricCurve::getLength(const unsigned int n) const {
	assert (n > 0);
	return getDiscretization(n)->getPolygonLength();
}

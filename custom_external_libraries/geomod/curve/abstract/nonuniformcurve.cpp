#include "nonuniformcurve.h"

using namespace geomodcore;
using namespace std;


NonUniformCurve::~NonUniformCurve() 
{
	knots.clear();
}

const Knots& NonUniformCurve::getKnotValuesConstRef() const 
{
	return knots;
}

Knots& NonUniformCurve::getKnotValuesRef() 
{
	return knots;
}

Knots NonUniformCurve::getKnotValues() const
{
	return knots;
}

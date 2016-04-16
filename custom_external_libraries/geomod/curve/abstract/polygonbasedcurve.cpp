#include "polygonbasedcurve.h"

using namespace std;
using namespace geomodcore;


PolygonBasedCurve::PolygonBasedCurve(const Polygon& poly)
	: polygon(poly)
{ 
}

PolygonBasedCurve::PolygonBasedCurve() 
{ 
}


PolygonBasedCurve::~PolygonBasedCurve() 
{
}

PolygonBasedCurve::PolygonBasedCurve(const PolygonBasedCurve& curve)
	: polygon(curve.getPolygonConstRef())
{
}

const Polygon& PolygonBasedCurve::getPolygonConstRef() const 
{
	return polygon;
}


Polygon& PolygonBasedCurve::getPolygonRef() 
{
	return polygon;
}

Polygon PolygonBasedCurve::getPolygon() const
{
	return polygon;
}

unsigned int PolygonBasedCurve::getNrOfControlPoints() const 
{
	return (unsigned int) polygon.size();
}

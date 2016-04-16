#pragma once
#include "../geomodcore.h"
#include "bsplinecurve.h"


namespace geomodcore {

	class NURBSCurve : public BSplineCurve
	{

	public:

		// Konstruktoren
		NURBSCurve(const Polygon& poly, const Knots& knotValues, const unsigned int K = 3);
		NURBSCurve(const PolygonBasedCurve& curve, const Knots& knotValues, const unsigned int K = 3);
		NURBSCurve(const BSplineCurve& curve);

		NURBSCurve(const Polygon& poly, const unsigned int K = 3);
		NURBSCurve(const PolygonBasedCurve& curve, const unsigned int K = 3);

		// Standard-Konstruktor
		NURBSCurve(void);
		// Copy-Konstruktor
		NURBSCurve(const NURBSCurve& curve);
		// Destruktor
		~NURBSCurve(void);
		// Assignment operator
		NURBSCurve& operator= (const NURBSCurve& curve);

		GPoint getPoint(double t) const;
	
	};

}
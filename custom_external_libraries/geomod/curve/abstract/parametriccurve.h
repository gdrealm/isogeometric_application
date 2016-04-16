#include "curve.h"

#pragma once

namespace geomodcore {

	class ParametricCurve : public Curve
	{

	public:

		/**
		* Gibt den dreidimensionalen Punkt am Kurvenparameter t aus.
		*/
		virtual GPoint getPoint(double t) const = 0;

		/**
		* Gibt den minimalen und maximalen Parameterwert der Kurve aus
		*/
		virtual MinMaxP1d minmaxT() const = 0;

		Polygon* getDiscretization(const unsigned int n) const;

		double getLength(const unsigned int n) const;
	};
}
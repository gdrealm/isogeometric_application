#pragma once
#include "../geomodcore.h"
#include "abstract/polygonbasedcurve.h"


namespace geomodcore {

	/**
	* Lagrange-Kurven über einem Kontrollpolygon
	* 
	* \author asche
	*/
	class LagrangeCurve : public PolygonBasedCurve {

	public:
		// Konstruktor
		LagrangeCurve(const Polygon& controlPts);
		LagrangeCurve(const PolygonBasedCurve& curve);

		// Standard-Konstruktor
		LagrangeCurve(void);
		// Copy-Konstruktor
		LagrangeCurve(const LagrangeCurve& curve);
		// Destruktor
		~LagrangeCurve(void);
		// Assignment operator
		LagrangeCurve& operator= (const LagrangeCurve& curve);

		GPoint getPoint(double t) const;


		MinMaxP1d minmaxT() const;


		/*
		* Berechnet das Langrange-Polygon mit Index i am Parameter t
		* \param i Index
		* \param t Kurvenparameter
		*/
		static double lagrange(const unsigned int i, const double t, const unsigned int K);

		
		/**
		* Bestimmung der abgeleiteten Kurve.
		*/
		LagrangeCurve* getDerivation(void) const;

		/**
		* Berechnet die Bogenlänge der Kurve aus der Ableitung.
		*/
		double getArcLength(void) const;

	protected:
		void update() ;
	};
}
#pragma once
#include "../geomodcore.h"
#include "abstract/polygonbasedcurve.h"


namespace geomodcore {

	class BezierCurve : public PolygonBasedCurve
	{

	public:

		// Konstruktor
		BezierCurve(const Polygon& controlPts);
		BezierCurve(const PolygonBasedCurve& curve);

		// Standard-Konstruktor
		BezierCurve(void);
		// Copy-Konstruktor
		BezierCurve(const BezierCurve& curve);
		// Destruktor
		~BezierCurve(void);
		// Assignment operator
		BezierCurve& operator= (const BezierCurve& curve);


	public:

		/**
		* Berechnung des Funktionswerts des i-ten Bernstein-Polynoms vom Grad n an der Stelle t
		* \param i Index
		* \param n Grad
		* \param t Stelle
		*/
		static double bernstein(const unsigned int i, const unsigned int n, const double t) {
			return double(boost::math::binomial_coefficient<double>(n, i) * pow(t, (double) i) * pow(1.0 - t, (double) (n - i))); 
		}

		/**
		* \brief Graderhöhung
		* Erhöht den Grad der Bezier-Kurve um 1 und fügt einen weiteren Kontrollpunkt ein
		*/
		void increaseDegree(void);

		GPoint getPoint(double t) const;

		// Gibt den minimalen und maximalen Parameterwert der Kurve aus
		MinMaxP1d minmaxT() const;
		
		/**
		* Bestimmung der abgeleiteten Kurve.
		*/
		BezierCurve* getDerivation(void) const;

		/**
		* Berechnet die Bogenlänge der Kurve aus der Ableitung.
		*/
		double getArcLength(void) const;
					
	};
}
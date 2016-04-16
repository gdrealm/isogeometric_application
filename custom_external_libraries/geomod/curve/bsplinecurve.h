#pragma once

#include "../geomodcore.h"
#include "abstract/polygonbasedcurve.h"
#include "abstract/nonuniformcurve.h"



namespace geomodcore {

	class BSplineCurve : public PolygonBasedCurve, public NonUniformCurve
	{

	protected:
		unsigned int K;

	public:

		/**
		* Konstruktor mit Grad K und optionalem Knotenvektor
		*/ 
		BSplineCurve(const Polygon& poly, const unsigned int K = 3);
		BSplineCurve(const PolygonBasedCurve& curve, const unsigned int K = 3);
		BSplineCurve(const Polygon& poly, const Knots& knotValues, const unsigned int K = 3);
		BSplineCurve(const PolygonBasedCurve& curve, const Knots& knotValues, const unsigned int K = 3);

		// Standard-Konstruktor
		BSplineCurve(void);
		// Copy-Konstruktor
		BSplineCurve(const BSplineCurve& curve);
		// Destruktor
		~BSplineCurve(void);
		// Assignment operator
		BSplineCurve& operator= (const BSplineCurve& curve);

		unsigned int getK() const;

		void setK(const unsigned int K);


		/**
		* Bestimmung der abgeleiteten Kurve.
		*/
		BSplineCurve* getDerivation(void) const;

		/**
		* Berechnet die Bogenlänge der Kurve aus der Ableitung.
		*/
		double getArcLength(void) const;


		GPoint getPoint(double t) const;


		// Gibt den minimalen und maximalen Parameterwert der Kurve aus
		MinMaxP1d minmaxT() const;


		/**
		* Berechnet rekursiv den Wert der B-Spline-Ansatzfunktion
		* \param u Kurvenparameter
		* \param i Index im Knotenvektor (Kontrollpunkt)
		* \param r Rekursionstiefe
		* \param K Grad des Polynoms
		* \param N Anzahl der Punkte
		* \param knots Knotenvektor der Länge N+K+1
		* \return Wert der Ansatzfunktion
		*/
		static double BSplineFunction(const double u, const unsigned int i, const unsigned int r, const unsigned int N, const unsigned int K, const Knots& knots);

		/**
		* Berechnet rekursiv den Wert der B-Spline-Ansatzfunktion bei uniformen Knotenvektor
		* \param u Kurvenparameter
		* \param i Index im Knotenvektor (Kontrollpunkt)
		* \param r Rekursionstiefe
		* \param K Grad des Polynoms
		* \param N Anzahl der Punkte
		* \return Wert der Ansatzfunktion
		*/
		static double BSplineFunctionUniform(const double u, const unsigned int i, const unsigned int r, const unsigned int N, const unsigned int K);

		/*
		* Initialisiert den Knotenvektor, so dass die K-1 äusseren Knoten mehrfach gesetzt sind, um die 
		* Endpunkt-Interpolation sicherzustellen.
		*/
		static Knots initUniformKnotValues(const unsigned int n, const unsigned int K);
	
	};
}
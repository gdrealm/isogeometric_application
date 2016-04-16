#pragma once
#include "../../geometry/polygon.h"

namespace geomodcore {

	class Curve 
	{

	public:
		/**
		* Erzeugt eine äquidistante Diskretisierung der Kurve.
		* Die Sampling-Größe darf nicht negativ sein. 
		* \param n Abtastschritte
		*/
		virtual Polygon* getDiscretization(const unsigned int n) const = 0;

		/**
		* Berechnet eine Näherung für die Länge der Kurve anhand einer Diskretisierung der Kurve
		* mit Parameterbereich L und der daraus resultierenden chordalen Länge.
		* Die Methode kann überschrieben werden, sofern eine vollständige Berechnung der Bogenlänge möglich ist.
		* \param n Anzahl der Schritte
		*/
		virtual double getLength(const unsigned int n) const = 0;

		/**
		* Bestimmung der abgeleiteten Kurve.
		*/
		Curve* getDerivation(void) const { };

		/**
		* Berechnet die Bogenlänge der Kurve aus der Ableitung.
		*/
		virtual double getArcLength(void) const = 0;
					
	};
}
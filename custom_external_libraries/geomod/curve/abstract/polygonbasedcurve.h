#pragma once

#include "../../geomodcore.h"
#include "parametriccurve.h"

namespace geomodcore {

	/**
	* \brief Polygon-basierte Kurven
	*
	* Virtuelle Klasse für Polygon-basierte Kurven
	*/
	class PolygonBasedCurve : public ParametricCurve
	{

	protected:
		// Polygon mit Kontrollpunkten
		Polygon polygon;

		// Destruktor
		~PolygonBasedCurve(void);

	public:
		// Konstruktor mit Polygonzug
		PolygonBasedCurve(const Polygon& poly);
		// Standard-Konstruktor
		PolygonBasedCurve(void);
		// Copy-Konstruktor
		PolygonBasedCurve(const PolygonBasedCurve& curve);
		// Assignment operator
		PolygonBasedCurve& operator= (const PolygonBasedCurve& curve);


		/**
		* \brief Konstantes Kontrollpolygon
		* Gibt das Kontrollpolygon der Kurve als konstante Referenz zurück
		* \return Kontrollpolygon 
		*/
		const Polygon& getPolygonConstRef(void) const;
			
		/**
		* \brief Kontrollpolygon
		* Gibt das Kontrollpolygon der Kurve als Referenz zurück
		* \return Kontrollpolygon 
		*/
		Polygon& getPolygonRef(void);

		/**
		* \brief Kontrollpolygon
		* Gibt das Kontrollpolygon der Kurve zurück
		* \return Kontrollpolygon 
		*/
		Polygon getPolygon(void) const;

		/*
		* \brief Anzahl der Kontrollpunkte
		* Gibt die Anzahl der Kontrollpunkte zurück (size_t als unsigned int)
		* \return Anzahl
		*/
		unsigned int getNrOfControlPoints(void) const;

	};

}
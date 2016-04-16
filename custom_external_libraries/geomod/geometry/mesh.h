#pragma once

#include "../geomodcore.h"
#include "polygon.h"

namespace geomodcore {

	class Mesh {

	public:

		/**
		* \brief Erzeugt Polygon-Liste zum Mesh
		* Transformiert das Kontrollnetz in eine Liste aus Polygonen. 
		* Hierbei wird jedes Face in ein einzelnes, im UZS laufendes Polygon gewandelt, so dass die
		* Spezialisierung der Klasse dann wahlweise Dreieck-, Vier- oder n-Eck ist.
		* \return Polygon-Array
		*/
		virtual std::vector<Polygon> toPolygons(void) const = 0;

		/**
		* \brief Löschen
		* Abstrakte Methode zum Löschen des Netzes.
		*/
		virtual void clear(void) = 0;

		/**
		* \brief minimale XYZ-Koordinaten
		* Bestimmt die minimalsten und maximalsten Koordinaten im Netz 
		* \return minimale XYZ-Koordinaten
		*/
		virtual MinMaxG minmaxG(void) const = 0;

		/**
		* \brief affine XYZ-Transformation
		* Ermöglicht es, eine homogene affine Transformation im Koordinatenbereich des Kontrollnetzes durchzuführen
		* \param T affine Transformationsmatrix
		*/
		virtual void transformG(const AffTransformG& T) = 0;

		/**
		* \brief minimale Bild-Koordinaten
		* \return minimale UV-Koordinaten
		*/
		virtual MinMaxI minmaxI(void) const = 0;


		/**
		* \brief affine UV-Transformation
		* Ermöglicht es, eine homogene affine Transformation im Bildbereich des Kontrollnetzes durchzuführen
		* \param T affine Transformationsmatrix
		*/
		virtual void transformI(const AffTransformP& T) = 0;
	};

}

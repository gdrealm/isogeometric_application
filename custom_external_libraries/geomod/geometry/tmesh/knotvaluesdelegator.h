#pragma once

#include "../../geomodcore.h"
#include "../tmesh.h"

namespace geomodcore {

	class TMesh;

	class KnotValuesDelegator
	{

	private:
		// TMesh
		const TMesh* const mesh;
		// Point-Location-Objekt
		const NaiveTMeshPointLocation& pl;
	public:

		KnotValuesDelegator(const TMesh* const mesh, const NaiveTMeshPointLocation& pl);

		
		/**
		* \brief Herleitung der Knotenvektoren aus dem TMesh
		* Bestimmt die lokalen horizontalen & vertikalen Knotenwerte aus dem aktuellen Zustand des TMesh am Vertex v.
		* \return horizontale & vertikale Knotenwerte
		*/
		KnotsN getKnotValues(const SVertex* const v) const;
		KnotsN getKnotValues(const Point2d& i) const;

	private:
		
		std::vector<Segment2d> getIntersectionSegments(Point2d) const;

		void trimKnotMaps(std::vector<std::map<double,double>>& maps) const;

	};

}
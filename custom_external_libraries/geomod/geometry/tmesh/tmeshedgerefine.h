#pragma once

#include "../../geomodcore.h"
#include "simpletmesh.h"

namespace geomodcore {

	typedef std::pair<TKnotsN,GPoint> BlendingFunctionEntry;
	typedef std::multimap<Vertex<>*, BlendingFunctionEntry> BlendingFunctionMap;

	// Alter Vertex, neuer Vertex und Blending-Factor
	typedef boost::tuple<Vertex<>*, Vertex<>*, const double> BFMatrixEntry;
	typedef std::list<BFMatrixEntry> BlendingFunctionMatrix;	

	class TMeshEdgeRefine
	{

	private:
		// SimpleTMesh
		SimpleTMesh<>* const tmesh;

		BlendingFunctionMap bfMap;
		BlendingFunctionMatrix bfm;

	public:

		TMeshEdgeRefine(SimpleTMesh<>* const tmesh);
		~TMeshEdgeRefine();

		/**
		* \brief Edge-Refinement
		* Fügt einen neuen Bildpunkt i mit Parameterkoordinaten p auf der übergebenen Kante ein
		* Es muss sichergestellt werden, dass der Parameterpunkt innerhalb des gleichen SimpleTMesh-Objekts wie i liegt!
		*/
		void refineEdge(Edge* oldEdge, const IPoint& iPt, const PPoint& pPt);

	private:
		Vertex<>* splitEdge(Edge* const oldEdge, const IPoint& iPt, const PPoint& pPt);

		bool handleViolation1();
		bool handleViolation2();
		
		/**
		* \brief Geometrie-Phase des Refining
		* Führt die Geometrie-Phase des Refining auf den BlendingFunctions durch
		*/
		void updateGeometry();

		void blendingRefine(const double k, Vertex<>* oldVertex, const unsigned int component);

		// Vergleich von Knotenwerte
		bool compareKnots(const TKnotsN& knotsA, const TKnotsN& knotsB) const;

		void printBfMap() const;
		void printBfMatrix() const;

	};

}
#pragma once

#include "../../geomodcore.h"
#include "simpletmesh.h"

namespace geomodcore {

	class TMeshFaceRefine
	{

	private:
		// SimpleTMesh
		SimpleTMesh<>* const tmesh;
		// Point-Location-Objekt
	public:

		TMeshFaceRefine(SimpleTMesh<>* const tmesh);


		/**
		* \brief Face-Refinement
		* Fügt einen neuen Bildpunkt i ein
		*/
		void refineFace(Face* face, const IPoint& iPt, const PPoint& pPt);

		Edge* splitFace(Face* oldF, Vertex<>* v0, Vertex<>* v1);
	};

}
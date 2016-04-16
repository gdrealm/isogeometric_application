#include "tmeshfacerefine.h"

using namespace std;
using namespace geomodcore;

TMeshFaceRefine::TMeshFaceRefine(SimpleTMesh<>* const tmesh) : tmesh(tmesh)
{
}

Edge* TMeshFaceRefine::splitFace(Face* oldF, Vertex<>* v0, Vertex<>* v1)
{
	assert(v0 != nullptr);
	assert(v1 != nullptr);
	assert(v0 != v1);

	// leere Faces
	Face* newF1 = new Face();
	Face* newF2 = new Face();

	// neue Kante
	Edge* newE = new Edge();
	newE->vertices[0] = v0;
	newE->vertices[1] = v1;

	newF1->edges.insert(newE);
	newF2->edges.insert(newE);

	// Bestimme Menge der Kanten der beiden neuen Faces
	bool horizontalSplit = (v0->p[1] == v1->p[1]) && (v0->p[0] != v1->p[0]);
	for (set<Edge*>::iterator eIt = oldF->edges.begin(); eIt != oldF->edges.end(); eIt++)
	{
		if (horizontalSplit) {
			if (max((*eIt)->vertices[0]->p[1], (*eIt)->vertices[1]->p[1]) < v0->p[1]) newF1->edges.insert(*eIt);
			else newF2->edges.insert(*eIt);
		}
		else {
			if (max((*eIt)->vertices[0]->p[0], (*eIt)->vertices[1]->p[0]) < v0->p[0]) newF2->edges.insert(*eIt);
			else newF2->edges.insert(*eIt);
		}
	}

	// Ersetze altes Face durch neue Faces
	//	faces.erase(oldF);
	//	faces.insert(newF1);
	//	faces.insert(newF2);

	return newE;
}


void TMeshFaceRefine::refineFace(Face* const face, const IPoint& iPt, const PPoint& pPt)
{
	// füge neuen Punkt auf Kante ein, wenn ein Strahl vom Punkt p die Kante in eine der vier Richtungen schneidet

	// Suche Punktkoordinaten für Refine
	const Vertex<>* v;
	
	const MinMaxI& minmaxI = SimpleTMesh<>::minmax(face).first;

	// Süd
	const IPoint southI = { iPt[0],minmaxI.first[1] };
	const IPoint southP = tmesh->parameterPoint(southI);
	// nur Einfügen, falls Bildpunkt nicht vorhanden ist
	v = tmesh->findVertexConst(southI);
	if (v != nullptr)
	{
		tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(southI, southP));
	}

	// Nord
	const IPoint northI = { iPt[0],minmaxI.second[1] };
	const PPoint northP = tmesh->parameterPoint(northI);
	v = tmesh->findVertexConst(northI);
	if (v != nullptr)
	{
		tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(northI, northP));
	}

	// Center
	tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(iPt, pPt));

	// West
	const IPoint westI = { minmaxI.first[0], iPt[1] };
	const PPoint westP = tmesh->parameterPoint(westI);
	v = tmesh->findVertexConst(westI);
	if (v != nullptr)
	{
		tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(westI, westP));
	}

	// Ost
	const IPoint eastI = { minmaxI.second[0], iPt[1] };
	const PPoint eastP = tmesh->parameterPoint(eastI);
	v = tmesh->findVertexConst(westI);
	if (v != nullptr)
	{
		tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(eastI, eastP));
	}
}
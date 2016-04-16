#include "tmeshedgerefine.h"

using namespace std;
using namespace geomodcore;


TMeshEdgeRefine::TMeshEdgeRefine(SimpleTMesh<>* const tmesh) : tmesh(tmesh)
{
}

TMeshEdgeRefine::~TMeshEdgeRefine()
{
	bfMap.clear();
	bfm.clear();
}

Vertex<>* TMeshEdgeRefine::splitEdge(Edge* const oldEdge, const IPoint& iPt, const PPoint& pPt)
{

	// Split-Vertex
	PPoint p = {(oldEdge->vertices[0]->p[0] + oldEdge->vertices[1]->p[0]) / 2, (oldEdge->vertices[0]->p[1] + oldEdge->vertices[1]->p[1]) / 2};
	IPoint i = {(oldEdge->vertices[0]->i[0] + oldEdge->vertices[1]->i[1]) / 2, (oldEdge->vertices[0]->i[1] + oldEdge->vertices[1]->i[1]) / 2};
	GPoint x = {0,0,0};
	Vertex<>* newV = new Vertex<>(p,i,x,1);

	// neue Kanten
	Edge* newE1 = new Edge();
	newE1->vertices[0] = oldEdge->vertices[0];
	newE1->vertices[1] = newV;

	Edge* newE2 = new Edge();
	newE2->vertices[0] = newV;
	newE1->vertices[1] = oldEdge->vertices[1];

	// Suche Faces mit alter Kante
	std::set<Face*> faces = tmesh->getFaces();
	for (set<Face*>::iterator fIt = faces.begin(); fIt != faces.end(); fIt++)
	{
		Face* f = *fIt;
		set<Edge*>::iterator foundEdge = f->edges.find(oldEdge);
		if (foundEdge != f->edges.end())
		{
			f->edges.erase(foundEdge);
			f->edges.insert(newE1);
			f->edges.insert(newE2);
		}
	}
	return newV;
}


void TMeshEdgeRefine::refineEdge(Edge* oldEdge, const IPoint& iPt, const PPoint& pPt)
{
	splitEdge(oldEdge,iPt, pPt);

	// Wir arbeiten auf blendingFunctions, d.h. Zeiger auf die Knotenwerte und die Abbildung auf die Vertices
	BlendingFunctionMap bfMap;
	BlendingFunctionMatrix bfm;

	// Einschränkung auf die benachbarten Vertices, da diese einzig von den Violations betroffen sind
	handleViolation1();
	// Violations 2 nur auf den veränderten BlendingFunctions prüfen
	handleViolation2();

	// Geometrie-Phase durchführen
	updateGeometry();
}

bool TMeshEdgeRefine::handleViolation1()
{
	bool violationFound = false;
	for (set<Vertex<>*>::iterator vIt = tmesh->getVertices().begin(); vIt != tmesh->getVertices().end(); vIt++)
	{
		Vertex<>* v = *vIt;

		// Mengendifferenz bestimmen
		vector<double> diffH;
		TKnotsN derivedKnots = tmesh->getKnotValues(v);
		const TKnotsN& bf = v->knots;
		if (bf.size() < 2) continue;				// Violation 3

		std::set_difference(derivedKnots[0].begin(), derivedKnots[0].end(), bf[0].begin(), bf[0].end(), inserter(diffH, diffH.begin()));
		assert(diffH.size() <= 1);
		// Falls die Violation am Rand des Knotenvektors erfolgt, dann ignoriere diese
		if ((diffH.size() == 1) && (diffH.front() > bf[0].front()) && (diffH.front() < bf[0].back())) 
		{
			cout << "Horizontale Violation 1 am Knoten (" << v->i[0] << "," << v->i[1] << ")" << endl;
			double k = diffH.front();
			blendingRefine(k, v, 0);
			violationFound = true;
		} 

		vector<double> diffV;
		std::set_difference(derivedKnots[1].begin(), derivedKnots[1].end(), bf[1].begin(), bf[1].end(), inserter(diffV, diffV.begin()));
		assert(diffV.size() <= 1);
		if ((diffV.size() == 1) && (diffV.front() > bf[1].front()) && (diffV.front() < bf[1].back())) 
		{
			cout << "Vertikale Violation 1 am Knoten (" << v->i[0] << "," << v->i[1] << ")" << endl;
			double k = diffV.front();
			blendingRefine(k, v, 1);
			violationFound = true;
		} 
	}
	cout << "Violation 1 - Prüfung abgeschlossen." << endl;
	return violationFound;
}


bool TMeshEdgeRefine::handleViolation2()
{
	std::set<Vertex<>*> uniqueVertices;
	//			printBfMap(bfMap);
	for (BlendingFunctionMap::const_iterator keyIt = bfMap.cbegin(); keyIt != bfMap.cend(); keyIt++)
	{
		// aktueller Vertex
		uniqueVertices.insert(keyIt->first);
		// Knotenwerte zurücksetzen
		keyIt->first->knots.fill(TKnots());
	}
	bool violationFound = false;
	for (set<Vertex<>*>::const_iterator vIt = uniqueVertices.cbegin(); vIt != uniqueVertices.cend(); vIt++)
	{
		Vertex<>* v = *vIt;
		pair<BlendingFunctionMap::const_iterator, BlendingFunctionMap::const_iterator> bfAtVertex;
		bfAtVertex = bfMap.equal_range(v);

		TKnotsN derivedKnots = tmesh->getKnotValues(v);
		// horizontale Violation 2
		for (BlendingFunctionMap::const_iterator bfIt = bfAtVertex.first; bfIt != bfAtVertex.second; bfIt++)
		{
			Vertex<>* v = bfIt->first;
			BlendingFunctionEntry bfEntry = bfIt->second;
			TKnotsN& bfKnots = bfEntry.first;

			vector<double> diffH;
			vector<double> diffV;

			std::set_difference(bfKnots[0].begin(), bfKnots[0].end(), derivedKnots[0].begin(), derivedKnots[0].end(), inserter(diffH, diffH.begin()));
			std::set_difference(bfKnots[1].begin(), bfKnots[1].end(), derivedKnots[1].begin(), derivedKnots[1].end(), inserter(diffV, diffV.begin()));

			// falls abweichende Knoten vorliegen (also Violation 2), dann erzeuge neue Kontrollpunkte

			if ((diffH.size() > 0) || (diffV.size() > 0)) {
				// Falls die Violation am Rand des Knotenvektors erfolgt, dann ignoriere diese
				if (diffH.size() > 0)
				{
					for (vector<double>::const_iterator kIt = diffH.cbegin(); kIt != diffH.cend(); kIt++)
					{
						double k = *kIt;
						//const double k = std::max(*kIt, 0.);
						if ((k > derivedKnots[0].front()) && (k < derivedKnots[0].back()))
						{
							cout << "Horizontale Violation 2 am Knoten " << "(" << v->i[0] << "," << v->i[1] << ")" << endl;
							PPoint newPt = { k, bfKnots[1][(tmesh->getK() + 1)/2] };
							v->knots = bfKnots;
							IPoint iPt = tmesh->imagePoint(newPt);		// TODO hier können ggf mehrdeutige Bildkoordinaten entstehen
							cout << "Neuer Kontrollpunkt bei i = ("<< iPt[0] << "," << iPt[1] << "), p=(" << newPt[0] << "," << newPt[1] << ")" << endl;
							tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(iPt, newPt));
							violationFound = true;

							continue;
						} else 
						{
							cout << "Ignorierte horizontale Violation 2 am Knoten " << "(" << v->i[0] << "," << v->i[1] << ")" << endl;
							continue;
						}
					}
				}
				else if (diffV.size() > 0)
				{
					for (vector<double>::const_iterator kIt = diffV.cbegin(); kIt != diffV.cend(); kIt++)
					{
						double k = *kIt;
						//const double k = std::max(*kIt, 0.);
						if ((k > derivedKnots[1].front()) && (k < derivedKnots[1].back()))
						{
							cout << "Vertikale Violation 2 am Knoten " << "(" << v->i[0] << "," << v->i[1] << ")" << endl;
							PPoint newPt = { bfKnots[0][(tmesh->getK() + 1)/2], k };
							v->knots = bfKnots;
							IPoint iPt = tmesh->imagePoint(newPt);		// TODO hier können ggf mehrdeutige Bildkoordinaten entstehen
							cout << "Neuer Kontrollpunkt bei i = ("<< iPt[0] << "," << iPt[1] << "), p=(" << newPt[0] << "," << newPt[1] << ")" << endl;
							tmesh->insertionPoints.push_back(pair<IPoint,PPoint>(iPt, newPt));
							violationFound = true;

							continue;
						} else 
						{
							cout << "Ignorierte vertikale Violation 2 am Knoten " << "(" << v->i[0] << "," << v->i[1] << ")" << endl;
							continue;
						}
					}
				}
			}	else {
				// Keine Violation 2, d.h. die Knotenwerte sind gesichert
				// diese werden nur bei leeren Knotenwerten überschrieben, da sonst die Violation-2-Knotenwerte überschrieben werden
				if (v->knots.size() == 0)
				{
					v->knots = bfKnots;
				}
			}
		}
	}
	cout << "Violation 2 - Prüfung abgeschlossen." << endl;
	return violationFound;
}


void TMeshEdgeRefine::blendingRefine(const double k, Vertex<>* oldVertex, const unsigned int component)
{
	const unsigned int K = tmesh->getK();
	// Wandeln auf STL-Vector (also KnotsN), da mehr als K+2 Elemente
	KnotsN oldBF = { Knots(oldVertex->knots[0].cbegin(), oldVertex->knots[0].cend()), Knots(oldVertex->knots[1].cbegin(), oldVertex->knots[1].cend()) }  ;		
	
	oldBF[component].push_back(k);
	sort(oldBF[component].begin(), oldBF[component].end());
	
	unsigned int posK = 0;
	for (Knots::const_iterator it = oldBF[component].cbegin(); it != oldBF[component].cend(); it++)
	{
		if (*it == k) break;
		posK++;
	}
	assert((oldBF[component].size() == K + 3) && (oldBF[1 - component].size() == 5));
	assert(posK >= 0 && posK < K + 3);

	// gesplittete Blending-Functions 
	TKnotsN splitBF1, splitBF2;
	splitBF1.fill(TKnots());
	splitBF2.fill(TKnots());
	for (int i = 0; i < K+2; i++)
	{
		splitBF1[1 - component][i] = oldBF[1 - component][i];
		splitBF2[1 - component][i] = oldBF[1 - component][i];
	}

	// Knotenwert einfügen und durchsortieren
	const Knots& knots = oldBF[component];

	// aufgesplittete BlendingFunctions
	Knots::const_iterator first = knots.cbegin();
	Knots::const_iterator last = knots.cbegin() + 5;
	for (int i = 0; i < K+2; i++)
	{
		splitBF1[component][i] = knots[i];
		splitBF2[component][i] = knots[i+1];

	}

	// Blending-Vorfaktoren
	const double c = (k - knots[0]) / (knots[4] - knots[0]);
	const double d = (knots[5] - k) / (knots[5] - knots[1]);

	// Suche gesplittete Blending-Functions anhand der Position von k im Knotenvektor, dies muss auch bei mehrfachen Knoten funktionieren
	Vertex<>* v1It;
	Vertex<>* v2It;
	if (posK <= (K+1)/2)  
	{
		v1It = tmesh->findNextVertex(oldVertex, component == 0 ? 3 : 2);
		v2It = oldVertex;
	} else {
		v1It = oldVertex;
		v2It = tmesh->findNextVertex(oldVertex, component == 0 ? 1 : 0);
	}
	assert(v1It != nullptr);
	assert(v2It != nullptr);
	assert(v1It != v2It);
	GPoint gPoint1 = { oldVertex->p[0] * c, oldVertex->p[1] * c, oldVertex->p[2] * c };
	GPoint gPoint2 = { oldVertex->p[0] * d, oldVertex->p[1] * d, oldVertex->p[2] * d };
	pair<TKnotsN,GPoint> pair1(splitBF1, GPoint(gPoint1));
	pair<TKnotsN,GPoint> pair2(splitBF2, GPoint(gPoint2));

	bfMap.insert(pair<Vertex<>*,pair<TKnotsN,GPoint>>(v1It, pair1));
	bfMap.insert(pair<Vertex<>*,pair<TKnotsN,GPoint>>(v2It, pair2));
	bfm.push_back(BFMatrixEntry(oldVertex, v1It, c));
	bfm.push_back(BFMatrixEntry(oldVertex, v2It, d));
}

void TMeshEdgeRefine::updateGeometry()
{
	// Abbildung vom Vertex auf die Geometrie mit BlendingFactor
	std::map<Vertex<>*, std::pair<GPoint, double>> blendedGeometry;
	for (BlendingFunctionMatrix::iterator entryIt = bfm.begin(); entryIt != bfm.end(); entryIt++) 
	{
		BFMatrixEntry& entry = *entryIt;
		Vertex<>* const vOld = get<0>(entry);
		Vertex<>* const vNew = get<1>(entry);
		const double factor = get<2>(*entryIt);

		// Suche nach neuem Vertex für Geometrie-Summierung
		std::map<Vertex<>*, std::pair<GPoint, double>>::iterator resultIt = blendedGeometry.find(vNew);
		if (resultIt != blendedGeometry.cend())
		{
			// Vertex schon vorhanden
			std::pair<GPoint, double>& pair = resultIt->second;
			GPoint gPt = { pair.first[0] + vOld->p[0] * factor, pair.first[1] + vOld->p[1] * factor, pair.first[2] + vOld->p[2] * factor };
			pair.first = GPoint(gPt);
			pair.second = pair.second + factor;
		} else 
		{
			// Vertex neu einfügen
			GPoint gPt = { vOld->p[0] * factor, vOld->p[1] * factor, vOld->p[2] * factor };
			std::pair<GPoint, double> pair(gPt, factor);
			blendedGeometry.insert(std::pair<Vertex<>* const, std::pair<GPoint, double>> (vNew, pair));
		}
	}

	// Skalierung um summierte Faktoren durchführen und Geometrie-Koordinaten zurückschreiben
	for (std::map<Vertex<>*, std::pair<GPoint, double>>::iterator geomIt = blendedGeometry.begin(); geomIt != blendedGeometry.end(); geomIt++) 
	{
		GPoint gPt = { geomIt->second.first[0] / geomIt->second.second, geomIt->second.first[1] / geomIt->second.second, geomIt->second.first[2] / geomIt->second.second };
		geomIt->first->x = gPt;
	}
	return;

	// Lösche die neu ermittelten Geometrie-Punkte...
	for (BlendingFunctionMap::const_iterator bfIt = bfMap.cbegin(); bfIt != bfMap.cend(); bfIt++)
	{
		Vertex<>* const v = const_cast<Vertex<>* const>(bfIt->first);
		GPoint gPt = { 0, 0, 0};
		v->x = gPt;
	}
	// ... und setze die gesammelten Informationen der BFMap auf die Vertices
	for (BlendingFunctionMap::const_iterator bfIt = bfMap.cbegin(); bfIt != bfMap.end(); bfIt++)
	{
		Vertex<>* const v = const_cast<Vertex<>* const>(bfIt->first);
		GPoint gPt = { v->p[0] + bfIt->second.second[0], v->p[1] + bfIt->second.second[1], v->p[2] + bfIt->second.second[2] };
		v->x = gPt;
	}
}

bool TMeshEdgeRefine::compareKnots(const TKnotsN& knotsA, const TKnotsN& knotsB) const
{
	assert(knotsA.size() == knotsB.size());
	for (unsigned int i = 0; i < knotsA.size(); i++)
	{
		if (equal(knotsA[i].cbegin(), knotsA[i].cend(), knotsB[i].cbegin()) == false) return false;
	}
	return true;
}

void TMeshEdgeRefine::printBfMap() const
{
	for (BlendingFunctionMap::const_iterator it = bfMap.cbegin(); it != bfMap.cend(); it++)
	{
		const Vertex<>* const v = it->first;
		cout << "VertexPtr " << v << " (" << v->i[0] << "," << v->i[1] << ")" << endl;
		pair<TKnotsN,GPoint> bf = it->second;
		cout << "Knoten [";
		for (TKnots::const_iterator kIt = bf.first[0].cbegin(); kIt != bf.first[0].cend(); kIt++) cout << *kIt << " ";
		cout << "]" << endl <<"       [";
		for (TKnots::const_iterator kIt = bf.first[1].cbegin(); kIt != bf.first[1].cend(); kIt++) cout << *kIt << " ";
		cout << "]" << endl;
		cout << "Geometrie " << bf.second[0] << " " << bf.second[1] << " " << bf.second[2] << " "<< endl << endl;
	}
	cout << endl;
}


void TMeshEdgeRefine::printBfMatrix() const
{
	for (BlendingFunctionMatrix::const_iterator it = bfm.cbegin(); it != bfm.cend(); it++)
	{
		BFMatrixEntry entry = *it;
		cout << "alter Vertex (" << get<0>(entry)->i[0] << "," << get<0>(entry)->i[1] << ")" << endl;
		cout << "neuer Vertex (" << get<1>(entry)->i[0] << "," << get<1>(entry)->i[1] << ")" << endl;
		cout << "BlendingFactor " << get<2>(entry) << endl;
	}
	cout << endl;
}

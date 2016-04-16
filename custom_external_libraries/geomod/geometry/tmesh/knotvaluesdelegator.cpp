#include "tmeshfacerefine.h"

using namespace std;
using namespace geomodcore;

KnotValuesDelegator::KnotValuesDelegator(const TMesh* const mesh, const NaiveTMeshPointLocation& pl) : mesh(mesh), pl(pl)
{
}

KnotsN KnotValuesDelegator::getKnotValues(const Point2d& i) const
{
	CGAL::Object obj = pl.locate(i);

	// Handle werden benötigt, um den gefundenen Location-Query-Objekt aufzunehmen
	TMesh::Vertex_const_handle vh;

	if (CGAL::assign(vh, obj)) {
		return getKnotValues(vh.ptr());
	} else return KnotsN();
}

KnotsN KnotValuesDelegator::getKnotValues(const TMesh::Vertex* const v) const
{
	const Point2d& p = v->point();

	// Horizontale Schnitte mit Edges bestimmen
	const int K = mesh->getK();
	const int k = (K + 1) / 2;

	// Schnitt-Segmente, die abzusuchenden Koordinaten liegen maximal (K+1)/2 Bildkoordinaten entfernt
	std::vector<Segment2d> intersectionSegments = getIntersectionSegments(p);

	// Schnitt-Ergebnisse
	std::vector<CGAL::Object> intersectionsObjs[4];

	// Gefundene Knotenwerte werden als auf die x- bzw y-Mittelwerts der Halbkante abgelegt. Dies verhindert doppelt gefundene Werte bei entgegengesetzten Halbkanten!
	std::vector<std::map<double,double>> knotMaps;

	// Ablauf: Übernehme Koordinaten des Startpunkts der gefundenen Halbkanten
	for (int i = 0; i < 4; i++)
	{
		knotMaps.push_back(std::map<double,double>());
		CGAL::zone(const_cast<TMesh&>(*mesh), intersectionSegments[i], std::back_inserter(intersectionsObjs[i]), pl);
		for (std::vector<CGAL::Object>::const_iterator objIt = intersectionsObjs[i].cbegin(); objIt != intersectionsObjs[i].cend(); objIt++)
		{
			CGAL::Object obj = *objIt;
			TMesh::Vertex_handle hv;
			TMesh::Halfedge_handle he;
			TMesh::Face_handle hf;
			if (CGAL::assign(he, *objIt))
			{
				knotMaps[i].insert(pair<double,double>(he->source()->point()[i / 2], he->source()->data().vertex.st[i / 2]));
				knotMaps[i].insert(pair<double,double>(he->target()->point()[i / 2], he->target()->data().vertex.st[i / 2]));
			} else if (CGAL::assign(hv, *objIt))
			{
				knotMaps[i].insert(pair<double,double>(hv->point()[i / 2], hv->data().vertex.st[i / 2]));
			} 
		}
	}

	// Trimmen der gefunden Knotenwerte
	trimKnotMaps(knotMaps);

	assert(knotMaps[0].size() == k);
	assert(knotMaps[1].size() == k+1);
	assert(knotMaps[2].size() == k);
	assert(knotMaps[3].size() == k+1);

	// Kopiere horizontale und vertikale Knotenwerte in finales Array
	KnotsN knots(2);
	for (map<double, double>::const_iterator it = knotMaps[0].cbegin(); it != knotMaps[0].cend(); it++) {
		knots[0].push_back(it->second);
	}
	for (map<double, double>::const_iterator it = knotMaps[1].cbegin(); it != knotMaps[1].cend(); it++) {
		knots[0].push_back(it->second);
	}
	for (map<double, double>::const_iterator it = knotMaps[2].cbegin(); it != knotMaps[2].cend(); it++) {
		knots[1].push_back(it->second);
	}
	for (map<double, double>::const_iterator it = knotMaps[3].cbegin(); it != knotMaps[3].cend(); it++) {
		knots[1].push_back(it->second);
	}

	assert(knots[0].size() == K+2);
	assert(knots[1].size() == K+2);

	return knots;
}


std::vector<Segment2d> KnotValuesDelegator::getIntersectionSegments(Point2d p) const
{
	const int K = mesh->getK();
	const int k = 2000 * ((K + 1) / 2);		// TODO ggf absetzen

	std::vector<Segment2d> segments;

	segments.push_back(Segment2d(Point2d(p[0] - k, p[1]), p));		// West
	segments.push_back(Segment2d(p, Point2d(p[0] + k, p[1])));		// Ost
	segments.push_back(Segment2d(Point2d(p[0], p[1] - k), p));		// Süd
	segments.push_back(Segment2d(p, Point2d(p[0], p[1] + k)));		// Nord
	return segments;
}


void KnotValuesDelegator::trimKnotMaps(std::vector<std::map<double,double>>& maps) const
{
	const unsigned int k = (mesh->getK() + 1) / 2;

	// Falls Länge < ((K+1)/2+1) = k+1  (also 3 für für K = 3), dann fülle nicht gefundene Knotenwerte auf
	// Sonst schneide Elemente der Map heraus, so dass Gesamtlänge 3 entsteht

	// West
	if (maps[0].size() < k+1)
	{
		while (maps[0].size() < k+1) 
		{
			maps[0].insert(pair<double,double>(maps[0].begin()->first - 1, maps[0].begin()->second));
		}
	} else if (maps[0].size() > k+1) {
		while (maps[0].size() > k+1) 
		{
			maps[0].erase(maps[0].begin());
		}
	}

	// letzten Eintrag löschen aus West-Map
	maps[0].erase(--maps[0].end());

	// Ost
	if (maps[1].size() < k+1)
	{
		while (maps[1].size() < k+1) 
		{
			maps[1].insert(pair<double,double>(maps[1].rbegin()->first + 1, maps[1].rbegin()->second));
		}
	} else if (maps[1].size() > k+1) {
		while (maps[1].size() > k+1) 
		{
			maps[1].erase(--maps[1].end());
		}
	}

	// Süd
	if (maps[2].size() < k+1)
	{
		while (maps[2].size() < k+1) 
		{
			maps[2].insert(pair<double,double>(maps[2].begin()->first - 1, maps[2].begin()->second));
		}
	} else if (maps[2].size() > k+1) {
		while (maps[2].size() > k+1) 
		{
			maps[2].erase(maps[2].begin());
		}
	}
	// letzten Eintrag löschen aus Süd-Map
	maps[2].erase(--maps[2].end());

	// Nord
	if (maps[3].size() < k+1)
	{
		while (maps[3].size() < k+1) 
		{
			maps[3].insert(pair<double,double>(maps[3].rbegin()->first + 1, maps[3].rbegin()->second));
		}
	} else if (maps[3].size() > k+1) {
		while (maps[3].size() > k+1) 
		{
			maps[3].erase(--maps[3].end());
		}
	}
}

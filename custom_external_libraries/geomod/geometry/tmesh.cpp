#include "tmesh.h"

using namespace std;
using namespace geomodcore;

int TMesh::count = 0;

TMesh::TMesh(void) 
	: K(3)
{
	pl.attach(*this);
}

TMesh::TMesh(const TMesh& tmesh) 
	: TMeshArrangement(tmesh), K(3)
{
	pl.attach(*this);
}


TMesh::~TMesh(void) 
{
	pl.detach();
}

TMesh& TMesh::operator=(const TMesh& tmesh)
{
	return *(new TMesh(tmesh));
}

TMesh::TMesh(const std::vector<PPoint3d*>& points, const unsigned int rows, const unsigned int cols) 
	: K(3)
{
	assert(rows >= 2);
	assert(cols >= 2);

	// für Skalierung der Bildkoordinaten im Interval [0,1]

	// falls der vector leer ist, dann mit leeren Punkten auffüllen
	if (points.empty()) {
		(const_cast<vector<PPoint3d*>& >(points)).assign(cols*rows, new PPoint3d());
	}
	cout << "Erzeuge Netz (1/3)";

	// Segmente erzeugen und in leeres TMesh einfügen
	for (unsigned int row = 0; row < rows; row++) {
		Halfedge_handle heH = this->insert_in_face_interior(Segment2d(PPoint2d(0,row), PPoint2d(1,row)), this->unbounded_face());
		for (unsigned int col = 1; col < cols - 1; col++) {
			heH = this->insert_from_left_vertex(Segment2d(PPoint2d(col, row), PPoint2d(col+1, row)), heH->target());
		}
		cout << ".";
	}
	cout << endl << "Erzeuge Netz (2/3)";
	list<Segment2d> segments;
	for (unsigned int col = 0; col < cols; col++) {
		for (unsigned int row = 0; row < rows - 1; row++) {
			CGAL::insert(*this, Segment2d(PPoint2d(col, row), PPoint2d(col,row+1)));
			//segments.push_back(Segment2d(PPoint2d(col,row), PPoint2d(col,row+1)));
		}
		cout << ".";
	}
	//CGAL::insert(*this, segments.cbegin(), segments.cend());

	// Vertex-Injections (Geometrie- und Parameterkoordinaten)
	cout << endl << "Erzeuge Netz (3/3)";
	for (TMeshArrangement::Vertex_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint2d& iPt = vIt->point();
		const int u = (int) iPt[0].doubleValue();		// TODO Fehlerquelle!!!
		const int v = (int) iPt[1].doubleValue();
		vIt->data().vertex.p = *points[v * cols + u];
		vIt->data().vertex.st = PPoint2d(iPt[0],iPt[1]);
		vIt->data().vertex.w = 1;
		vIt->data().vertex.knots.assign(2, Knots());

		// Knotenwerte einfügen
		const int kD = (K+1) / 2;
		for (int k = -kD; k <= kD; k++) {
			double d = min(cols-1,(unsigned int) max(0,u+k));
			double e = min(rows-1,(unsigned int) max(0,v+k));
			vIt->data().vertex.knots[0].push_back(d);
			vIt->data().vertex.knots[1].push_back(e);
		}

		assert (vIt->data().vertex.knots[0].size() == K + 2);
		assert (vIt->data().vertex.knots[1].size() == K + 2);
	}
	cout << " fertig" << endl;

	// ID-Injection
	for (TMeshArrangement::Edge_iterator eIt = this->edges_begin(); eIt != this->edges_end(); eIt++) {
		//		eIt->data().edge.id = NEXT_EDGE_ID;
		eIt->twin()->data().edge.id = NEXT_EDGE_ID++;
	}
	pl.attach(*this);
}

unsigned int TMesh::getK() const
{
	return K;
};

void TMesh::setK(unsigned int K)
{
	assert (K % 2 == 1);
	this->K = K;
}

TMesh* TMesh::fromFile(const char* const filename)
{
	cout << "Importiere TMesh aus " << filename << endl;

	ifstream file;
	file.open(filename);
	string line;
	std::vector<std::string> linesplit;

	do 
	{
		getline(file, line);
	} while (line.empty());

	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));

	const int K = boost::lexical_cast<int>(linesplit[1]);
	cout << "> Grad K = " << K << endl;
	TMesh* tmesh = new TMesh();
	tmesh->setK(K);

	do 	{
		getline(file, line);
	} while (line.empty());

	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	const int nrOfVertices = boost::lexical_cast<int>(linesplit[1]);
	cout << "> Anzahl Vertices = " << nrOfVertices << endl;

	// Vertex-IDs auf Vertices abbilden 
	std::map<int, Vertex_handle> vertices;
	do 	{
		getline(file, line);
	} while (line.empty());

	while (line.empty() == false)
	{
		boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
		assert((int) linesplit.size() >= (9 + 2 * (K + 2)));		// id, u, v, s, t, x, y, z, w, 2x Knoten

		PPoint2d iPt = PPoint2d(boost::lexical_cast<double>(linesplit[1]), boost::lexical_cast<double>(linesplit[2]));
		Vertex_handle v = CGAL::insert_point(*tmesh, iPt, tmesh->pl);

		v->data().vertex.id = boost::lexical_cast<int>(linesplit[0]);
		v->data().vertex.st = PPoint2d(boost::lexical_cast<double>(linesplit[3]), boost::lexical_cast<double>(linesplit[4]));
		v->data().vertex.p = PPoint3d(boost::lexical_cast<double>(linesplit[5]), boost::lexical_cast<double>(linesplit[6]), boost::lexical_cast<double>(linesplit[7]));
		v->data().vertex.w = boost::lexical_cast<double>(linesplit[8]);

		const int nrOfKnotValues = (linesplit.size() - 9) / 2;

		Knots knotsH, knotsV;
		for (int i = 0; i < nrOfKnotValues; i++) 
		{
			knotsH.push_back (boost::lexical_cast<double>(linesplit[9 + i]));
			knotsV.push_back (boost::lexical_cast<double>(linesplit[9 + nrOfKnotValues + i]));
		}
		v->data().vertex.knots.push_back(knotsH);
		v->data().vertex.knots.push_back(knotsV);
		vertices.insert(pair<int,Vertex_handle>(v->data().vertex.id, v));
		getline(file, line);
	}

	// Segmente
	do {
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	const int nrOfEdges = boost::lexical_cast<int>(linesplit[1]);
	cout << "> Anzahl Kanten = " << nrOfEdges << endl;

	do {
		getline(file, line);
	} while (line.empty());
	while (line.empty() == false)
	{
		boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
		assert(linesplit.size() == 4);		// id, v0, v1, d
		int id0 = boost::lexical_cast<int>(linesplit[1]);
		int id1 = boost::lexical_cast<int>(linesplit[2]);
		Vertex_handle v0 = vertices.find(id0)->second;
		Vertex_handle v1 = vertices.find(id1)->second;

		Segment2d segment(v0->point(), v1->point());
		CGAL::insert(*tmesh, segment, tmesh->pl);
		getline(file, line);
	}

	file.close();

	if (tmesh->checkKnotValueConsistency() == false)
	{
		cout << "Achtung! Die importierten Knotenwerte sind falsch!" << endl << "Es werden neue Knotenwerte berechnet." << endl;
		for (TMeshArrangement::Vertex_iterator vIt = tmesh->vertices_begin(); vIt != tmesh->vertices_end(); vIt++)
		{
			vIt->data().vertex.knots = tmesh->getKnotValues(vIt.ptr());
		}
		assert(tmesh->checkKnotValueConsistency());
	}
	return tmesh;
}


void TMesh::reassignIds(void) 
{
	NEXT_VERTEX_ID = 0;
	NEXT_EDGE_ID = 0;
	NEXT_FACE_ID = 0;
	for (TMeshArrangement::Vertex_iterator it = this->vertices_begin(); it != this->vertices_end(); it++) {
		it->data().vertex.id = NEXT_VERTEX_ID++;
	}
	for (TMeshArrangement::Halfedge_iterator it = this->halfedges_begin(); it != this->halfedges_end(); it++) {
		it->data().edge.id = NEXT_EDGE_ID++;
	}
	for (TMeshArrangement::Face_iterator it = this->faces_begin(); it != this->faces_end(); it++) {
		it->data().face.id = NEXT_FACE_ID++;
	}
}

void TMesh::toFile(const char* const filename) const {

	const_cast<TMesh*>(this)->reassignIds();
	//if (this->checkKnotValueConsistency() == false) cout << "Achtung! Es wird ein inkonsistentes Netz gespeichert." << endl;
	cout << "Speichere TMesh als " << filename << endl;
	ofstream file;
	file.open (filename);

	// Grad des TMesh
	file << "K " << this->K << endl << endl;

	// vertices
	file << "vertices " << this->number_of_vertices() << endl;

	for (TMeshArrangement::Vertex_const_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint2d& iPt = vIt->point();
		const TVertex& vertex = vIt->data().vertex;

		file << vertex.id << " " << CGAL::to_double(iPt[0]) << " " << CGAL::to_double(iPt[1]) << " " << vertex.st[0] << " " << vertex.st[1] << " " << vertex.p[0]  << " " << vertex.p[1] << " " << vertex.p[2] << " " << vertex.w << " ";

		// Sind Knotenwerte vorhanden?
		if (vertex.knots.size() != 2) continue;

		for (Knots::const_iterator dIt = vertex.knots[0].cbegin(); dIt != vertex.knots[0].cend(); dIt++) {
			file << *dIt << " ";
		}
		for (Knots::const_iterator dIt = vertex.knots[1].cbegin(); dIt != vertex.knots[1].cend(); dIt++) {
			file << *dIt << " ";
		}
		file << endl;
	}

	// edges
	file << endl << "edges " << this->number_of_edges() << endl;
	for (TMeshArrangement::Edge_const_iterator eIt = this->edges_begin(); eIt != this->edges_end(); eIt++) {
		const int id = eIt->data().edge.id;
		const double d_squared = CGAL::squared_distance(eIt->source()->data().vertex.st, eIt->target()->data().vertex.st);
		const double d = sqrt(d_squared).doubleValue();
//		const double d = sqrt(CGAL::to_double(d_squared.exact()));

		const int v0 = eIt->source()->data().vertex.id;
		const int v1 = eIt->target()->data().vertex.id;

		file << id << " " << v0 << " " << v1 << " " << d << endl; 
	}

	// faces
	file << endl << "faces " << (this->number_of_faces() - 1) << endl;
	for (TMeshArrangement::Face_const_iterator fIt = this->faces_begin(); fIt != this->faces_end(); fIt++) {
		// Auslassen des unbounded-Randface
		if (fIt->is_unbounded() == false) {
			const int id = fIt->data().face.id - 1;		// -1 als Korrektur für das äußere unbound-Face

			file << id << " "; 
			// Circulator
			TMeshArrangement::Ccb_halfedge_const_circulator circ = fIt->outer_ccb();
			TMeshArrangement::Ccb_halfedge_const_circulator curr = circ;
			do {
				// Wähle minimale Halbkanten-ID als Edge-ID
				const int halfEdgeID = min(curr->data().edge.id, curr->twin()->data().edge.id);
				file << halfEdgeID << " ";
			} while (++curr != circ);
			file << endl;
		}
	}
	file << endl;
	file.close();
}

MinMaxI2d TMesh::minmaxUV(void) const
{
	double minu = MAXDOUBLE, minv = MAXDOUBLE, maxu = MINDOUBLE, maxv = MINDOUBLE;
	for (TMeshArrangement::Vertex_const_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint2d& iPt = vIt->point();
		minu = min(minu, iPt[0]);
		minv = min(minv, iPt[0]);
		maxu = max(maxu, iPt[1]);
		maxv = max(maxv, iPt[1]);
	}
	return MinMaxI2d(PPoint2d(minu,minv),PPoint2d(maxu,maxv));
}

MinMaxP2d TMesh::minmaxST(void) const
{
	double mins = MAXDOUBLE, mint = MAXDOUBLE, maxs = MINDOUBLE, maxt = MINDOUBLE;
	for (TMeshArrangement::Vertex_const_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint2d& p = vIt->data().vertex.st;
		mins = min(mins, p[0]);
		mint = min(mint, p[0]);
		maxs = max(maxs, p[1]);
		maxt = max(maxt, p[1]);
	}
	return MinMaxP2d(PPoint2d(mins,mint),PPoint2d(maxs,maxt));
}

void TMesh::transformUV(const AffTransform2d& T)
{
	for (TMesh::Vertex_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++)
	{
		Vertex_handle vh = this->modify_vertex(vIt, vIt->point());
		vh->point() = T.transform(vh->point());
	}
	return;


	// 1. Erzeuge leeres TMesh
	TMesh tmpTMesh;
	// 2. durchlaufe Edges des alten Mesh:
	// - transformiere Start- und Endpunkte
	// - füge diese in tmp ein
	list<Segment2d> segments;
	for (TMesh::Edge_const_iterator eIt = this->edges_begin(); eIt != this->edges_end(); eIt++) {
		PPoint2d p(eIt->source()->point());
		T.transform(p);
		PPoint2d q(eIt->target()->point());
		T.transform(q);
		segments.push_back(Segment2d(p,q));
	}
	CGAL::insert(tmpTMesh, segments.begin(), segments.end());

	// 3. leere dieses (this) T-Mesh und fülle es mit transformiertem Netz
	this->clear();
	this->assign(tmpTMesh);
	tmpTMesh.clear();
}

std::vector<Polygon>* TMesh::toPolygons(void) const
{
	std::vector<Polygon>* polygons = new std::vector<Polygon>;
	for (TMesh::Face_const_iterator fIt = this->faces_begin(); fIt != this->faces_end(); fIt++) {
		Polygon p;
		if (fIt->is_unbounded() == false) {

			// Circulator
			TMeshArrangement::Ccb_halfedge_const_circulator circ = fIt->outer_ccb();
			// erster Anfangspunkt der half-edges
			p.push_back(new PPoint3d(circ->source()->data().vertex.p));
			TMeshArrangement::Ccb_halfedge_const_circulator curr = circ;
			do {
				// fügt jeweils den Endpunkt der halfedge zum Polygon hinzu
				p.push_back(new PPoint3d(curr->target()->data().vertex.p));
			} while (++curr != circ);
			polygons->push_back(p);
		}
	}
	return polygons;
}

MinMaxXYZ TMesh::minmaxXYZ(void) const 
{
	double minx = MAXDOUBLE, miny = MAXDOUBLE, minz = MAXDOUBLE, maxx = MINDOUBLE, maxy = MINDOUBLE, maxz = MINDOUBLE;
	for (TMeshArrangement::Vertex_const_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint3d& p = vIt->data().vertex.p;
		minx = min(minx, p.hx());
		miny = min(miny, p.hy());
		minz = min(minz, p.hz());
		maxx = max(maxx, p.hx());
		maxy = max(maxy, p.hy());
		maxz = max(maxz, p.hz());
	}
	return MinMaxXYZ(PPoint3d(minx,miny,minz),PPoint3d(maxx,maxy,maxz));
}

void TMesh::transformXYZ(const AffTransform3d& T)
{
	for (TMeshArrangement::Vertex_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		vIt->data().vertex.p = T.transform(vIt->data().vertex.p);
	}
}

void TMesh::transformST(const AffTransform2d& T)
{
	for (TMeshArrangement::Vertex_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		vIt->data().vertex.st = T.transform(vIt->data().vertex.st);
	}
}

void TMesh::clear(void)
{
	((TMeshArrangement*) this)->clear();
}

TMesh* TMesh::getSampleTMesh()
{
	return TMeshSamples::getSampleTMesh();
}

TMesh* TMesh::getSampleTMesh2()
{
	return TMeshSamples::getSampleTMesh2();
}

TMesh* TMesh::getSampleTMeshMoscow()
{
	return TMeshSamples::getSampleTMeshMoscow();
}


KnotsN TMesh::getKnotValues(const TMesh::Vertex* const v) const
{
	KnotValuesDelegator delegator(const_cast<const TMesh* const>(this), pl);
	return delegator.getKnotValues(v);
}

KnotsN TMesh::getKnotValues(const PPoint2d i) const
{
	KnotValuesDelegator delegator(const_cast<const TMesh* const>(this), pl);
	return delegator.getKnotValues(i);
}

void TMesh::updateKnotValues(void)
{
	for (TMesh::Vertex_iterator vIt = vertices_begin(); vIt != vertices_end(); vIt++) 
	{
		vIt->data().vertex.knots = getKnotValues(vIt.ptr());
	}
} 


std::pair<MinMaxI2d,MinMaxP2d> TMesh::minmax(const TMeshArrangement::Halfedge* const edge)
{
	const double& u1 = edge->source()->point()[0];
	const double& u2 = edge->target()->point()[0];
	const double& v1 = edge->source()->point()[1];
	const double& v2 = edge->target()->point()[1];
	const double& s1 = edge->source()->data().vertex.st[0];
	const double& s2 = edge->target()->data().vertex.st[0];
	const double& t1 = edge->source()->data().vertex.st[1];
	const double& t2 = edge->target()->data().vertex.st[1];
	return pair<MinMaxI2d, MinMaxP2d>(MinMaxI2d(PPoint2d(min(u1,u2),min(v1,v2)),PPoint2d(max(u1,u2),max(v1,v2))), MinMaxP2d(PPoint2d(min(s1,s2),min(t1,t2)),PPoint2d(max(s1,s2),max(t1,t2))));
}

std::pair<MinMaxI2d,MinMaxP2d> TMesh::minmax(const TMeshArrangement::Face* const face) 
{
	assert(face->is_unbounded() == false);
	double mins = MAXDOUBLE, mint = MAXDOUBLE, maxs = MINDOUBLE, maxt = MINDOUBLE;
	double minu = MAXDOUBLE, minv = MAXDOUBLE, maxu = MINDOUBLE, maxv = MINDOUBLE;

	TMeshArrangement::Ccb_halfedge_const_circulator curr = face->outer_ccb();
	TMeshArrangement::Ccb_halfedge_const_circulator start = curr;
	do {
		mins = min(mins, curr->source()->data().vertex.st[0]);
		mint = min(mint, curr->source()->data().vertex.st[1]);
		maxs = max(maxs, curr->source()->data().vertex.st[0]);
		maxt = max(maxt, curr->source()->data().vertex.st[1]);
		minu = min(minu, curr->source()->point()[0]);
		minv = min(minv, curr->source()->point()[1]);
		maxu = max(maxu, curr->source()->point()[0]);
		maxv = max(maxv, curr->source()->point()[1]);
	} while (++curr != start);
	return pair<MinMaxI2d, MinMaxP2d>(MinMaxI2d(PPoint2d(minu,minv),PPoint2d(maxu,maxv)), MinMaxP2d(PPoint2d(mins,mint),PPoint2d(maxs,maxt)));
}

void TMesh::refineByParameter(PPoint2d pPt) 
{
	set<PPoint2d> imagePoints = imagePoint(pPt);
	return refine(pPt, *imagePoints.cbegin());


	// Speichere zu jedem Face die Bildkoordinaten des Parameterpunkts
	std::map<TMeshArrangement::Face*, PPoint2d> refineFaces;
	// Suche Faces für Refine-Operation
	for (TMeshArrangement::Face_iterator faceIt = this->faces_begin(); faceIt != faces_end(); faceIt++)
	{
		if (faceIt->is_unbounded()) continue;
		const pair<MinMaxI2d,MinMaxP2d> minmax = TMesh::minmax(faceIt.ptr());
		const double& mins = minmax.second.first[0], mint = minmax.second.first[1], maxs = minmax.second.second[0], maxt = minmax.second.second[1];
		const double& minu = minmax.first.first[0], minv = minmax.first.first[1], maxu = minmax.first.second[0], maxv = minmax.first.second[1];

		if (mins <= pPt[0] && pPt[0] <= maxs && mint <= pPt[1] && pPt[1] <= maxt) 
		{ 
			const double u = minu + (pPt[0]-mins)/(maxs-mins) * (maxu-minu);
			const double v = minv + (pPt[1]-mint)/(maxt-mint) * (maxv-minv);
			refineFaces.insert(std::pair<TMeshArrangement::Face*, PPoint2d>(faceIt.ptr(), PPoint2d(u,v)));
		}
	}
	if (refineFaces.size() == 0) 
	{
		cout << "Es wurden keine Flächen gefunden, in denen der Parameterpunkt (" << pPt << ") liegt." << endl;
		return;
	} else if (refineFaces.size() == 1) 
	{
		return refine(pPt, refineFaces.cbegin()->second);
	} else 
	{
		cout << "Achtung! Die Bildkoordinaten zum Parameterpunkt (" << pPt << ") sind nicht eindeutig" << endl;
		// TODO
		return;
	}
}

void TMesh::refineByImage(PPoint2d iPt) 
{
	const PPoint2d pPt = parameterPoint(iPt);
	if (pPt.x() == pPt.x() && pPt.y() == pPt.y())
	{
		refine(parameterPoint(iPt), iPt);
	} else {
		cout << "Achtung! Die Parameterkoordinaten zum Bildpunkt (" << iPt << ") sind fehlerhaft." << endl;
	}
}

void TMesh::refine(PPoint2d pPt, PPoint2d iPt) 
{
	insertionPoints.push_back(pair<PPoint2d, PPoint2d>(iPt, pPt));
	refine();
}

void TMesh::refine(void) 
{
	while (!insertionPoints.empty())
	{
		cout << "Verbleibende Anzahl einzufügender Punkte: " << insertionPoints.size() << endl;
		const std::pair<PPoint2d, PPoint2d>& point = insertionPoints.front();
		const PPoint2d& iPt = point.first;
		const PPoint2d& pPt = point.second;

		std::set<TMeshArrangement::Vertex*> undefinedVertices;
		cout << "Starte Refinement an i = (" << iPt[0] << " " << iPt[1] << ") und p = (" << pPt << ")" << endl;

		// Handle werden benötigt, um den gefundenen Location-Query-Objekt aufzunehmen und das Refinement durchzuführen
		Vertex_const_handle vertexHandle;
		Halfedge_const_handle edgeHandle;
		Face_const_handle faceHandle;
		CGAL::Object obj = pl.locate(iPt);

		// Fallunterscheidung, ob i auf einer Kante liegt
		if (CGAL::assign(vertexHandle, obj)) {
			cout << "> Handle-Typ ist Vertex" << endl;
			cout << "** Kein Refinement nötig" << endl;
		} else if (CGAL::assign(edgeHandle, obj)) {
			cout << "> Handle-Typ ist Edge" << endl;
			// Refinement auf Edge ausführen
			refineEdge(const_cast<Halfedge* const>(edgeHandle.ptr()), iPt, pPt);
			insertMissingEdges(iPt, pPt);
			updateKnotValues();
		} else if (CGAL::assign(faceHandle, obj)) {
			// neue Punkte auf umliegenden Edges erzeugen und dort Refinement ausführen
			cout << "> Handle-Typ ist Face" << endl;
			refineFace(const_cast<Face* const>(faceHandle.ptr()), iPt, pPt);
		} else 
		{
			cout << "Achtung! Es wurde kein CGAL-Objekt zum Bildpunkt (" << iPt << ") lokalisiert.";
		}
		insertionPoints.pop_front();
	}
}

void TMesh::refineEdge(TMeshArrangement::Halfedge* const edge, const PPoint2d iPt, const PPoint2d pPt)
{
	TMeshEdgeRefine refine(this, pl);
	refine.refineEdge(edge, iPt, pPt);
}

void TMesh::refineFace(TMeshArrangement::Face* const face, const PPoint2d iPt, const PPoint2d pPt)
{
	TMeshFaceRefine refine(this, pl);
	refine.refineFace(face, iPt, pPt);
}


bool TMesh::checkKnotValueConsistency(void) const
{
	for (TMeshArrangement::Vertex_const_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++)
	{
		const KnotsN& localKnots = vIt->data().vertex.knots;
		const KnotsN derKnots = getKnotValues(vIt.ptr());

		if ((std::equal(localKnots[0].cbegin(), localKnots[0].cend(), derKnots[0].cbegin()) && std::equal(localKnots[1].cbegin(), localKnots[1].cend(), derKnots[1].cbegin())) == false)
		{
			return false;
		}
	}
	return true;
}

void TMesh::insertMissingEdges(const PPoint2d iPt, const PPoint2d pPt)
{
	// Liste der neuen Kante
	list<Segment2d> segments;

	// 1. Suche (Vertex-)Objekt zum Bildpunkt i
	CGAL::Object obj = pl.locate(iPt);
	Vertex_const_handle vertexHandle;
	Halfedge_const_handle edgeHandle;
	Halfedge_const_handle oppEdgeHandle;

	if (CGAL::assign(vertexHandle, obj) == false)
	{
		cout << "Achtung! Es wurde kein CGAL-Objekt zum Bildpunkt (" << iPt << ") lokalisiert.";
		return;
	}

	// 2. Bestimme anliegende Faces am Bildpunkt i
	set<Face* const> faces;
	TMeshArrangement::Halfedge_around_vertex_const_circulator heAtVertexStart = vertexHandle->incident_halfedges();
	TMeshArrangement::Halfedge_around_vertex_const_circulator heAtVertexCurr = heAtVertexStart;
	do {
		if (heAtVertexCurr->face()->is_unbounded() == false) 
		{
			faces.insert(const_cast<Face* const>(heAtVertexCurr->face().ptr()));
		}
	} while (heAtVertexStart != ++heAtVertexCurr);

	// 3. durchlaufe Faces zum Bildpunkt
	for (set<TMeshArrangement::Face* const>::iterator fIt = faces.begin(); fIt != faces.end(); fIt++) {
		MinMaxP2d minmax = TMesh::minmax(*fIt).second;
		const double& mins = minmax.first[0], maxs = minmax.second[0], mint = minmax.first[1], maxt = minmax.second[1];
		// 4. bestimme alle Vertices zum aktuellen Face
		set<const Vertex* const> vertices;
		TMeshArrangement::Ccb_halfedge_const_circulator start = (*fIt)->outer_ccb();	
		TMeshArrangement::Ccb_halfedge_const_circulator curr = start;
		do {
			TMeshArrangement::Ccb_halfedge_const_circulator second = curr;
			vertices.insert(curr->source().ptr());
		} while (++curr != start);	

		// 5. durchlaufe Vertices vom aktuellen Face
		for (set<const Vertex* const>::const_iterator faceVertexIt = vertices.cbegin(); faceVertexIt != vertices.cend(); faceVertexIt++)
		{
			const double& s = (*faceVertexIt)->data().vertex.st[0], t = (*faceVertexIt)->data().vertex.st[1];
			// 5.1 falls (s1 > mins) && (s1 < maxs) && (s1 == s2) && (t1 != t2) nehme neues Segment (v1,v2) auf
			if (  (s > mins)&&(s < maxs) && (s == pPt[0]) && (t != pPt[1]) )
			{
				segments.push_back(Segment2d(iPt, (*faceVertexIt)->point()));
			}
			// 5.2 horizontal analog
			if (  (t > mint)&&(t < maxt) && (t == pPt[1]) && (s != pPt[0]) )
			{
				segments.push_back(Segment2d(iPt, (*faceVertexIt)->point()));
			}
		}
	}

	// 6. neue Segment einfügen
	if (segments.size() > 0) 
	{
		CGAL::insert(*this, segments.begin(), segments.end());
	}
}

void TMesh::insertMissingEdges(void)
{
	// Durchlaufe alle Faces und vergleiche die Punkte paarweise
	vector<Segment2d> newSegments;

	// Bilde Mengen der Punkte (inkl. der isolated vertices!)
	for (TMeshArrangement::Face_iterator fIt = this->faces_begin(); fIt != this->faces_end(); fIt++)
	{
		if (fIt->is_unbounded()) continue;
		set<const Vertex* const> vertices;
		for (TMeshArrangement::Isolated_vertex_iterator it = fIt->isolated_vertices_begin(); it != fIt->isolated_vertices_end(); it++)
		{
			vertices.insert(it.ptr());
		}
		TMeshArrangement::Ccb_halfedge_const_circulator start = fIt->outer_ccb();	
		TMeshArrangement::Ccb_halfedge_const_circulator curr = start;
		do {
			TMeshArrangement::Ccb_halfedge_const_circulator second = curr;
			vertices.insert(curr->source().ptr());
		} while (++curr != start);	

		// Paarweise vergleichen
		for (set<const Vertex* const>::const_iterator v1It = vertices.begin(); v1It != vertices.end(); v1It++)
		{
			const Vertex* const v1 = *v1It;
			for (set<const Vertex* const>::const_iterator v2It = vertices.cbegin(); v2It != vertices.cend(); v2It++)
			{
				if (v1It == v2It) continue;
				const Vertex* const v2 = *v2It;
				// horizontal bzw vertikal
				if (((v1->data().vertex.st[1] == v2->data().vertex.st[1]) && (v1->data().vertex.st[0] != v2->data().vertex.st[0])) 
					|| ((v1->data().vertex.st[0] == v2->data().vertex.st[0]) && (v1->data().vertex.st[1] != v2->data().vertex.st[1])))
				{
					newSegments.push_back(Segment2d(v1->point(), v2->point()));
				}
			}
		}
	}
	CGAL::insert(*this, newSegments.begin(), newSegments.end());
}

TMeshArrangement::Halfedge* TMesh::findEdge(const TMeshArrangement::Vertex* const v, const Direction direction) const
{
	if (v == nullptr) return nullptr;
	TMeshArrangement::Halfedge_around_vertex_const_circulator heStart;
	TMeshArrangement::Halfedge_around_vertex_const_circulator heCurr;
	heStart = heCurr = v->incident_halfedges();
	do
	{
		TMeshArrangement::Halfedge* const he = const_cast<TMeshArrangement::Halfedge* const>(heCurr.ptr());
		switch(direction) {
		case NORTH:	
			{
				// Suche Source-Knoten als gegenüberliegendes Ende der inzidenten Halbkante
				if ( (heCurr->source()->point()[0] == v->point()[0]) && (heCurr->source()->point()[1] > v->point()[1])) return he;
				break;
			}
		case SOUTH:	
			{
				if ( (heCurr->source()->point()[0] == v->point()[0]) && (heCurr->source()->point()[1] < v->point()[1])) return he;
				break;
			}
		case EAST:	
			{
				if ( (heCurr->source()->point()[1] == v->point()[1]) && (heCurr->source()->point()[0] > v->point()[0])) return he;
				break;
			}
		case WEST:	
			{
				if ( (heCurr->source()->point()[1] == v->point()[1]) && (heCurr->source()->point()[0] < v->point()[0])) return he;
				break;
			}
		}
	} while (heStart != ++heCurr);
	return nullptr;
}

TMeshArrangement::Vertex* TMesh::findVertex(const TMeshArrangement::Vertex* const v, const Direction direction) const
{
	if (v == nullptr) return nullptr;
	Halfedge* hePtr = findEdge(v, direction);
	if (hePtr != nullptr)
	{
		return hePtr->source().ptr();
	} else return nullptr;
}

std::set<TMeshArrangement::Vertex*> TMesh::findNeighbourVertices(const TMeshArrangement::Vertex* const v, const unsigned int k) const
{
	set<TMeshArrangement::Vertex*> out;
	// Erzeuge Kopie des Pointers, da der Vertex später Rückgabe ist
	out.insert(const_cast<TMeshArrangement::Vertex*>(v));
	TMeshArrangement::Vertex* westV = const_cast<TMeshArrangement::Vertex*>(v);		
	TMeshArrangement::Vertex* eastV = const_cast<TMeshArrangement::Vertex*>(v);
	TMeshArrangement::Vertex* northV = const_cast<TMeshArrangement::Vertex*>(v);
	TMeshArrangement::Vertex* southV = const_cast<TMeshArrangement::Vertex*>(v);

	for (unsigned int i = 1; i <= k; i++)
	{
		westV = findVertex(westV, WEST);
		if (westV != nullptr) out.insert(westV);
		eastV = findVertex(eastV, EAST);
		if (eastV != nullptr) out.insert(eastV);
		northV = findVertex(northV, NORTH);
		if (northV != nullptr) out.insert(northV);
		southV = findVertex(southV, SOUTH);
		if (southV != nullptr) out.insert(southV);

	}
	return out;
}

std::set<PPoint2d> TMesh::imagePoint(const PPoint2d& p) const
{
	std::set<PPoint2d> points;
	for (TMeshArrangement::Face_const_iterator fIt = faces_begin(); fIt != faces_end(); fIt++)
	{
		if (fIt->is_unbounded()) continue;
		pair<MinMaxI2d, MinMaxP2d> minmax = TMesh::minmax(fIt.ptr());

		const double& i0 = minmax.first.first[0], i1 = minmax.first.second[0], j0 = minmax.first.first[1], j1 = minmax.first.second[1];
		const double& s0 = minmax.second.first[0], s1 = minmax.second.second[0], t0 = minmax.second.first[1], t1 = minmax.second.second[1];

		if ((s0 <= p[0]) && (p[0] <= s1) && (t0 <= p[1]) && (p[1] <= t1))
		{
			const double u = i0 + (i1-i0) * (p[0] - s0) / (s1-s0);
			const double v = j0 + (j1-j0) * (p[1] - t0) / (t1-t0);
			points.insert(PPoint2d(u,v));
		}
	}
	return points;
}


PPoint2d TMesh::parameterPoint(const PPoint2d& iPt) const
{
	CGAL::Object obj = pl.locate(iPt);
	Vertex_const_handle vertexHandle;
	Halfedge_const_handle edgeHandle;
	Face_const_handle faceHandle;
	PPoint2d pPt;
	if (CGAL::assign(vertexHandle, obj)) {
		return vertexHandle->data().vertex.st;
	} else if (CGAL::assign(edgeHandle, obj)) {
		const pair<MinMaxI2d,MinMaxP2d> minmax = TMesh::minmax(edgeHandle.ptr());
		const double& mins = minmax.second.first[0], mint = minmax.second.first[1], maxs = minmax.second.second[0], maxt = minmax.second.second[1];
		const double& minu = minmax.first.first[0], minv = minmax.first.first[1], maxu = minmax.first.second[0], maxv = minmax.first.second[1];
		// vertikale Kante?
		const double s = (minu == maxu ? mins : mins + (iPt[0] - minu)/(maxu-minu) * (maxs-mins));
		// horizontale Kante?
		const double t = (minv == maxv ? mint : mint + (iPt[1] - minv)/(maxv-minv) * (maxt-mint));
		pPt = PPoint2d(s,t);
	} else if (CGAL::assign(faceHandle, obj)) {
		if (faceHandle->is_unbounded()) 
		{
			cout << "Achtung! Der Bildpunkt (" << iPt << ") liegt außerhalb des TMesh.";
			pPt = PPoint2d(NAN,NAN);
		} else {
			const pair<MinMaxI2d,MinMaxP2d> minmax = TMesh::minmax(faceHandle.ptr());
			const double& mins = minmax.second.first[0], mint = minmax.second.first[1], maxs = minmax.second.second[0], maxt = minmax.second.second[1];
			const double& minu = minmax.first.first[0], minv = minmax.first.first[1], maxu = minmax.first.second[0], maxv = minmax.first.second[1];
			// vertikale Kante?
			const double s = (minu == maxu ? mins : mins + (iPt[0] - minu)/(maxu-minu) * (maxs-mins));
			// horizontale Kante?
			const double t = (minv == maxv ? mint : mint + (iPt[1] - minv)/(maxv-minv) * (maxt-mint));
			pPt = PPoint2d(s,t);
		}
	} else 
	{
		cout << "Achtung! Es wurde kein CGAL-Objekt zum Bildpunkt (" << iPt << ") lokalisiert.";
		pPt = PPoint2d(NAN,NAN);
	}
	cout << pPt << endl;
	return pPt;
}

void TMesh::scaleToSimplex()
{
	MinMaxXYZ minmax = minmaxXYZ();
	const double t0 = -minmax.first[0];
	const double t1 = -minmax.first[1];
	const double t2 = -minmax.first[2];
	const double s0 = minmax.second[0] - minmax.first[0] < EPSILON ? 1. : 1. / (minmax.second[0] - minmax.first[0]);
	const double s1 = minmax.second[1] - minmax.first[1] < EPSILON ? 1. : 1. / (minmax.second[1] - minmax.first[1]);
	const double s2 = minmax.second[2] - minmax.first[2] < EPSILON ? 1. : 1. / (minmax.second[2] - minmax.first[2]);
	AffTransform3d T(s0, 0, 0, s0*t0, 0, s1, 0, s1*t1, 0, 0, s2, s2*t2);
	for (TMeshArrangement::Vertex_iterator vIt = this->vertices_begin(); vIt != this->vertices_end(); vIt++) {
		const PPoint3d& p = vIt->data().vertex.p;
		vIt->data().vertex.p = PPoint3d(p[0] * s0 + t0 * s0, p[1] * s1 + t1 * s1, p[2] * s2 + t2 * s2);
	}
}

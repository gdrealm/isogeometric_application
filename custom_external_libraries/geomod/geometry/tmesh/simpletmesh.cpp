#include "simpletmesh.h"

using namespace std;
using namespace geomodcore;

template <unsigned int D, unsigned int K>
SimpleTMesh<>::SimpleTMesh(void) : m_K(K)
{
}

SimpleTMesh<>::SimpleTMesh(std::set<boost::shared_ptr<Vertex<> > > vertices) : m_K(K)
{
}

SimpleTMesh<>::SimpleTMesh(const SimpleTMesh& tmesh) : m_K(K)
{
    this->m_vertices.insert(tmesh.m_vertices.begin(), tmesh.m_vertices.end());
}


SimpleTMesh<>::~SimpleTMesh(void) 
{
}

SimpleTMesh<>& SimpleTMesh<>::operator=(const SimpleTMesh<>& tmesh)
{
    this->m_vertices.insert(tmesh.m_vertices.begin(), tmesh.m_vertices.end());
    return *this;
}

template <unsigned int D, unsigned int K> SimpleTMesh<D, K>::SimpleTMesh(const std::vector<GPoint>& points, boost::array<unsigned int, D> size) 
{
    for (unsigned int i = 0; i < D; i++)
        assert(size[i] >= 2);
    std::cout << "TODO" << std::endl;
}

unsigned int SimpleTMesh<>::getK() const
{
    return K;
};


std::set<boost::shared_ptr<Vertex<> > > SimpleTMesh<>::getVertices() 
{
    /*    std::set<Vertex<>*> vertices;
    for (set<Face*>::const_iterator fIt = faces.cbegin(); fIt != faces.cend(); fIt++)
    {
    for (set<Edge*>::const_iterator eIt = (*fIt)->edges.cbegin(); eIt != (*fIt)->edges.cend(); eIt++)
    {
    vertices.insert((*eIt)->vertices.cbegin(),    (*eIt)->vertices.cend());
    }
    }*/
    return m_vertices;
}


std::set<const boost::shared_ptr<Vertex<> > > SimpleTMesh<>::getVerticesConst() const
{
    std::set<const boost::shared_ptr<Vertex<> > > verticesConst;
    for (std::set<boost::shared_ptr<Vertex<> > >::const_iterator vIt = m_vertices.cbegin(); vIt != m_vertices.cend(); ++vIt) 
    {
        verticesConst.insert(*vIt);
    }
    return verticesConst;
}

std::set<boost::shared_ptr<Edge> > SimpleTMesh<>::getEdges()
{
    std::set<boost::shared_ptr<Edge> > edges;
    for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = m_vertices.begin(); vIt != m_vertices.end(); ++vIt)
    {
        edges.insert((*vIt)->edges.begin(), (*vIt)->edges.end());
    }
    return edges;
}

std::set<const boost::shared_ptr<Edge> > SimpleTMesh<>::getEdgesConst() const
{
    std::set<const boost::shared_ptr<Edge> > edges;
    for (std::set<boost::shared_ptr<Vertex<> > >::const_iterator vIt = m_vertices.cbegin(); vIt != m_vertices.cend(); ++vIt) 
    {
        edges.insert((*vIt)->edges.cbegin(), (*vIt)->edges.cend());
    }
    return edges;
}

std::set<const boost::shared_ptr<Face> > SimpleTMesh<>::getFacesConst() const
{
    std::set<const boost::shared_ptr<Face> > faces;
    std::set<const boost::shared_ptr<Edge> > edges;
    for (std::set<const boost::shared_ptr<Edge> >::const_iterator eIt = edges.cbegin(); eIt != edges.cend(); ++eIt)
    {
        faces.insert((*eIt)->faces.cbegin(), (*eIt)->faces.cend());
    }
    return faces;
}

template<> std::set<boost::shared_ptr<Face> > SimpleTMesh<>::getFaces()
{
    std::set<boost::shared_ptr<Face> > faces;
    std::set<boost::shared_ptr<Edge> > edges;
    for (std::set<boost::shared_ptr<Edge> >::iterator eIt = edges.begin(); eIt != edges.end(); ++eIt)
    {
        faces.insert((*eIt)->faces.cbegin(), (*eIt)->faces.cend());
    }
    return faces;
}

boost::shared_ptr<SimpleTMesh<> > SimpleTMesh<>::fromFile(const char* filename)
{
    cout << "Importiere TMesh aus " << filename << endl;

    // ID-Maps für importierte Daten, dies sind nicht die IDs in den neu erzeugten TMeshs
    std::map<int, boost::shared_ptr<Vertex<> > > vertices;
    std::map<int, boost::shared_ptr<Edge> > edges;

    std::ifstream file;
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

    do     {
        getline(file, line);
    } while (line.empty());

    boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
    const int nrOfVertices = boost::lexical_cast<int>(linesplit[1]);
    cout << "> Anzahl Vertices = " << nrOfVertices << endl;

    // Vertex-IDs auf Vertices abbilden 
    do     {
        getline(file, line);
    } while (line.empty());

    while (line.empty() == false)
    {
        boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
        assert((int) linesplit.size() >= (9 + 2 * (K + 2)));        // id, u, v, s, t, x, y, z, w, 2x Knoten

        int id = boost::lexical_cast<double>(linesplit[0]);
        IPoint iPt = { boost::lexical_cast<double>(linesplit[1]), boost::lexical_cast<double>(linesplit[2]) };
        PPoint pPt = { boost::lexical_cast<double>(linesplit[3]), boost::lexical_cast<double>(linesplit[4]) };
        GPoint xPt = { boost::lexical_cast<double>(linesplit[5]), boost::lexical_cast<double>(linesplit[6]), boost::lexical_cast<double>(linesplit[7]) };
        const double w = boost::lexical_cast<double>(linesplit[8]);

        boost::shared_ptr<Vertex<> > v = boost::shared_ptr<Vertex<> >(new Vertex<>(pPt, iPt, xPt, w));

        size_t nrOfKnotValues = (linesplit.size() - 9) / 2;

        for (int i = 0; i < nrOfKnotValues; i++) 
        {
            v->knots[0][i] = boost::lexical_cast<double>(linesplit[9 + i]);
            v->knots[1][i] = boost::lexical_cast<double>(linesplit[9 + nrOfKnotValues + i]);
        }

        vertices.insert(std::pair<int, boost::shared_ptr<Vertex<> > >(id, v));
        getline(file, line);
    }

    // Kanten
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
        assert(linesplit.size() == 4);        // id, v0, v1, d
        int id = boost::lexical_cast<int>(linesplit[0]);
        int vertex0Id = boost::lexical_cast<int>(linesplit[1]);
        int vertex1Id = boost::lexical_cast<int>(linesplit[2]);
        boost::shared_ptr<Vertex<> > v0 = vertices.find(vertex0Id)->second;
        boost::shared_ptr<Vertex<> > v1 = vertices.find(vertex1Id)->second;
        assert(v0 != 0);
        assert(v1 != 0);
        boost::shared_ptr<Edge> e = boost::shared_ptr<Edge>(new Edge());
        e->vertices[0] = v0;
        e->vertices[1] = v1;
        edges.insert(pair<int, boost::shared_ptr<Edge> >(id, e));

        // Einfügen der NSOW-Kante am Vertex
        if ((v0->i[0] < v1->i[0]) || ((v0->i[0] == v1->i[0]) && (v0->i[1] < v1->i[1])))
        {
            // Kantenorientierung v0->v1, bestimme Richtung der Kante
            if (v0->i[0] == v1->i[0])
            {
                // vertikal, also bei v0 Nordkante, bei v1 Südkante
                v0->edges[0] = e;
                v1->edges[2] = e;
            } else {
                // horizontal, also bei v0 Ostkante, bei v1 Westkante
                v0->edges[1] = e;
                v1->edges[3] = e;
            }
        } else {
            // Kantenorientierung v1->v0, bestimme Richtung der Kante
            if (v0->i[0] == v1->i[0])
            {
                // vertikal, also bei v1 Nordkante, bei v0 Südkante
                v0->edges[2] = e;
                v1->edges[0] = e;
            } else {
                // horizontal, also bei v1 Ostkante, bei v2 Westkante
                v0->edges[3] = e;
                v1->edges[1] = e;
            }
        }
        getline(file, line);
    }
    do {
        getline(file, line);
    } while (line.empty());

    // Faces
    while (line.empty() == false)
    {
        boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
        assert(linesplit.size() >= 5);    // 4 Edges + ID
        boost::shared_ptr<Face> f = boost::shared_ptr<Face>(new Face());
        int id = boost::lexical_cast<int>(linesplit[0]);    // unnötig
        for (int i = 1; i < linesplit.size(); i++)
        {
            int edgeId = boost::lexical_cast<int>(linesplit[i]);    
            boost::shared_ptr<Edge> edge = edges.find(edgeId)->second;
            f->edges.insert(edge);
            edge->faces.insert(f);
        }
        getline(file, line);
    }
    file.close();

    // T-Mesh erzeugen
    boost::shared_ptr<SimpleTMesh<> > tmesh = boost::shared_ptr<SimpleTMesh<> >(new SimpleTMesh<2,3>());
    for (std::map<int, boost::shared_ptr<Vertex<> > >::iterator vIt = vertices.begin(); vIt != vertices.end(); ++vIt)
    {
        tmesh->m_vertices.insert(vIt->second);
    }

    if (tmesh->checkKnotValueConsistency() == false)
    {
        cout << "Achtung! Die importierten Knotenwerte sind falsch!" << endl << "Es werden neue Knotenwerte berechnet." << endl;
        for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = tmesh->m_vertices.begin(); vIt != tmesh->m_vertices.end(); ++vIt)
        {
            TKnotsN knots = tmesh->getKnotValues(*vIt);
            (*vIt)->knots = knots;//.assign(knots);
            (*vIt)->knots[1] = knots[1];
        }
        assert(tmesh->checkKnotValueConsistency());
    }
    return tmesh;
}

void SimpleTMesh<>::toFile(const char* filename) const {}

void SimpleTMesh<>::updateKnotValues(void) {}

void SimpleTMesh<>::clear() {
    m_vertices.clear();
}

MinMaxG SimpleTMesh<>::minmaxG(void) const
{
    GPoint min, max;
    std::set<const boost::shared_ptr<Vertex<> > > vertices = this->getVerticesConst();
    for (std::set<const boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        const GPoint& x = (*vIt)->x;
        min[0] = std::min(min[0], x[0]);
        min[1] = std::min(min[1], x[1]);
        min[2] = std::min(min[2], x[2]);
        max[0] = std::max(max[0], x[0]);
        max[1] = std::max(max[1], x[1]);
        max[2] = std::max(max[2], x[2]);
    }
    return MinMaxG(min,max);
}

void SimpleTMesh<>::transformG(const AffTransformG& T)
{
    std::set<boost::shared_ptr<Vertex<> > > vertices = getVertices();
    for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = vertices.begin(); vIt != vertices.end(); ++vIt)
    {
        (*vIt)->x = transform((*vIt)->x,T);
    }
}

template <unsigned int D, unsigned int K>
std::pair<MinMaxI,MinMaxP> SimpleTMesh<D, K>::minmax(const boost::shared_ptr<Edge> edge)
{
    IPoint mini, maxi;
    PPoint minp, maxp;
    for (int i = 0; i < D; i++) {
        mini[i] = min(edge->vertices[0]->i[i],edge->vertices[1]->i[i]);
        maxi[i] = max(edge->vertices[0]->i[i],edge->vertices[1]->i[i]);
        minp[i] = min(edge->vertices[0]->p[i],edge->vertices[1]->p[i]);
        maxp[i] = max(edge->vertices[0]->p[i],edge->vertices[1]->p[i]);
    }
    return std::pair<MinMaxI,MinMaxP>(MinMaxI(mini,maxi), MinMaxP(minp,maxp));
}

template <unsigned int D, unsigned int K>
std::pair<MinMaxI,MinMaxP> SimpleTMesh<D, K>::minmax(const boost::shared_ptr<Face> face)
{
    IPoint mini, maxi;
    PPoint minp, maxp;
    mini.fill(MAXDOUBLE);
    maxi.fill(MINDOUBLE);
    minp.fill(MAXDOUBLE);
    maxp.fill(MINDOUBLE);

    std::set<boost::shared_ptr<Vertex<> > > vertices;
    for (std::set<boost::shared_ptr<Edge> >::const_iterator eIt = face->edges.cbegin(); eIt != face->edges.cend(); ++eIt)
    {
        vertices.insert((**eIt).vertices.cbegin(), (**eIt).vertices.cend());
    }

    for (std::set<boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        for (int i = 0; i < D; i++) {
            const Vertex<>& v = **vIt;
            mini[i] = min(mini[i], v.i[i]);
            maxi[i] = max(maxi[i], v.i[i]);
            minp[i] = min(minp[i], v.p[i]);
            maxp[i] = max(maxp[i], v.p[i]);
        }
    }
    return std::pair<MinMaxI,MinMaxP>(MinMaxI(mini,maxi), MinMaxP(minp,maxp));
}

template <unsigned int D, unsigned int K>
MinMaxI SimpleTMesh<D, K>::minmaxI(void) const
{
    IPoint mini, maxi;
    set<const boost::shared_ptr<Vertex<> > > vertices = this->getVerticesConst();
    for (set<const boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        for (int i = 0; i < D; i++) {
            mini[i] = min(mini[i], (*vIt)->i[i]);
            maxi[i] = max(maxi[i], (*vIt)->i[i]);
        }
    }
    return MinMaxI(mini,maxi);
}

void SimpleTMesh<>::transformI(const AffTransformI& T)
{
    std::set<boost::shared_ptr<Vertex<> > > vertices = getVertices();
    for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = vertices.begin(); vIt != vertices.end(); ++vIt)
    {
        (*vIt)->i = transform((*vIt)->i,T);
    }
}

template <unsigned int D, unsigned int K>
MinMaxP SimpleTMesh<D, K>::minmaxP(void) const
{
    PPoint minp, maxp;
    std::set<const boost::shared_ptr<Vertex<> > > vertices = this->getVerticesConst();
    for (std::set<const boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        for (int i = 0; i < D; i++) {
            minp[i] = min(minp[i], (*vIt)->p[i]);
            maxp[i] = max(maxp[i], (*vIt)->p[i]);
        }
    }
    return MinMaxP(minp,maxp);
}

void SimpleTMesh<>::transformP(const AffTransformP& T)
{
    std::set<boost::shared_ptr<Vertex<> > > vertices = getVertices();
    for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = vertices.begin(); vIt != vertices.end(); ++vIt)
    {
        (*vIt)->p = transform((*vIt)->p,T);
    }
}

template<> TKnotsN SimpleTMesh<>::getKnotValues(const boost::shared_ptr<Vertex<> > v) const
{
    cout << "TODO" << endl;
    return TKnotsN();
}


IPoint SimpleTMesh<>::imagePoint(const PPoint& pPt) const
{
    std::set<const boost::shared_ptr<Face> > faces = getFacesConst();
    for (std::set<const boost::shared_ptr<Face> >::const_iterator fIt = faces.cbegin(); fIt != faces.cend(); ++fIt)
    {
        std::pair<MinMaxI,MinMaxP> minmax = SimpleTMesh::minmax(*fIt);

        const double& i0 = minmax.first.first[0];
        const double& i1 = minmax.first.second[0];
        const double& j0 = minmax.first.first[1];
        const double& j1 = minmax.first.second[1];

        const double& s0 = minmax.second.first[0];
        const double& s1 = minmax.second.second[0];
        const double& t0 = minmax.second.first[1];
        const double& t1 = minmax.second.second[1];

        if ((s0 <= pPt[0]) && (pPt[0] <= s1) && (t0 <= pPt[1]) && (pPt[1] <= t1))        // TODO nD
        {
            IPoint i = {s0 + (s1-s0) * (pPt[0] - i0) / (i1-i0), t0 + (t1-t0) * (pPt[1] - j0) / (j1-j0)};
            return i;
        }
    }
    return IPoint();
}


PPoint SimpleTMesh<>::parameterPoint(const IPoint& iPt) const
{
    std::set<const boost::shared_ptr<Face> > faces = getFacesConst();
    for (std::set<const boost::shared_ptr<Face> >::const_iterator fIt = faces.cbegin(); fIt != faces.cend(); ++fIt)
    {
        std::pair<MinMaxI, MinMaxP> minmax = SimpleTMesh::minmax(*fIt);

        const double& i0 = minmax.first.first[0];
        const double& i1 = minmax.first.second[0];
        const double& j0 = minmax.first.first[1];
        const double& j1 = minmax.first.second[1];

        const double& s0 = minmax.second.first[0];
        const double& s1 = minmax.second.second[0];
        const double& t0 = minmax.second.first[1];
        const double& t1 = minmax.second.second[1];

        if ((i0 <= iPt[0]) && (iPt[0] <= i1) && (j0 <= iPt[1]) && (iPt[1] <= j1))        // TODO nD
        {
            PPoint p = {i0 + (i1-i0) * (iPt[0] - s0) / (s1-s0), j0 + (j1-j0) * (iPt[1] - t0) / (t1-t0)};
            return p;
            //points.insert(i);
        }
    }
    return PPoint();
}


std::vector<Polygon> SimpleTMesh<>::toPolygons(void) const
{
    // Erzeuge vereinfachte Polygon-Menge als Kanten zwischen den Punkten
    std::vector<Polygon> polygons;
    std::set<const boost::shared_ptr<Edge> > edges = getEdgesConst();
    for (std::set<const boost::shared_ptr<Edge> >::const_iterator eIt = edges.cbegin(); eIt != edges.cend(); ++eIt) {
        Polygon p;
        p.push_back((*eIt)->vertices[0]->x);
        p.push_back((*eIt)->vertices[1]->x);
        polygons->push_back(p);
    }
    return polygons;
}

const boost::shared_ptr<Vertex<> > SimpleTMesh<>::findVertexConst(const IPoint& iPt) const
{
    std::set<const boost::shared_ptr<Vertex<> > > vertices = this->getVerticesConst();
    for (std::set<const boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        if ((**vIt).i == iPt) return *vIt;
    }
    return nullptr;
}

boost::shared_ptr<Vertex<> > SimpleTMesh<>::findVertex(const IPoint& iPt)
{
    std::set<const boost::shared_ptr<Vertex<> > > vertices = this->getVerticesConst();
    for (std::set<const boost::shared_ptr<Vertex<> > >::const_iterator vIt = vertices.cbegin(); vIt != vertices.cend(); ++vIt)
    {
        if ((**vIt).i == iPt) return const_cast<Vertex<>*>(*vIt);
    }
    return nullptr;
}


template<> bool SimpleTMesh<>::checkKnotValueConsistency() const 
{ 
    for (std::set<boost::shared_ptr<Vertex<> > >::iterator vIt = m_vertices.begin(); vIt != m_vertices.end(); ++vIt)
    {
        const TKnotsN& localKnots = (*vIt)->knots;
        const TKnotsN derivedKnots = getKnotValues(*vIt);

        if (localKnots != derivedKnots) return false;
        /*if ((std::equal(localKnots[0].cbegin(), localKnots[0].cend(), derKnots[0].cbegin()) && std::equal(localKnots[1].cbegin(), localKnots[1].cend(), derKnots[1].cbegin())) == false)
        {
            return false;
        }*/
    }
    return true;

}


boost::shared_ptr<SimpleTMesh<> > SimpleTMesh<>::getSampleTMesh()
{
    const Grid& sampleGrid = Grid::getSampleGrid();
    return boost::shared_ptr<SimpleTMesh<> >(new SimpleTMesh<>(sampleGrid));
};


boost::shared_ptr<SimpleTMesh<> > SimpleTMesh<>::getSampleTMesh2(std::string filename)
{
    return SimpleTMesh<>::fromFile(filename.c_str());
};


boost::shared_ptr<Vertex<> > SimpleTMesh<>::findNextVertex(const boost::shared_ptr<Vertex<> > v, const unsigned int direction) const
{
    if (v == 0) return 0;
    boost::shared_ptr<Edge> e = findOutgoingEdge(v, direction);
    if (e != 0)
    {
        if (e->vertices[0] == v)
            return e->vertices[1];
        else
            return e->vertices[0];
    } else
        return 0;
}


template<> boost::shared_ptr<Edge> SimpleTMesh<>::findOutgoingEdge(const boost::shared_ptr<Vertex<> > v, const unsigned int inDirection) const
{
    if (v == 0) return 0;
    switch (inDirection) 
    {
    case 0 : // aus NORD nach SÜD
        return v->edges[2];
    case 1 : // aus WEST nach OST
        return v->edges[3];
    case 2 : // aus SÜD nach NORD
        return v->edges[0];
    case 3 : // aus OST nach WEST
        return v->edges[1];
    default:
        return nullptr;
    }
}

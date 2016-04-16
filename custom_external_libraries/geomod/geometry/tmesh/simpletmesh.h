#pragma once

#include "../../geomodcore.h"
#include "../grid.h"

namespace geomodcore {

    class Edge;
    class Face;

    template <unsigned int D = 2, unsigned int K = TSplineDegree>
    class Vertex 
    {
    private:
        static int STATIC_ID;
    public:
        boost::array<boost::shared_ptr<Edge>, (1 << D) > edges;
        GPoint x;
        PPoint p;
        PPoint i;
        double w;
        TKnotsN knots;
        const int ID;
        Vertex(const PPoint& p, const PPoint& i, const GPoint& x, const double w) : ID(STATIC_ID++)
        {
            this->p = p;
            this->i = i;
            this->x = x;
            this->w = w;
        }
    };

    class Edge
    {
    private:
        static int STATIC_ID;
    public:
        boost::array<boost::shared_ptr<Vertex<> >, 2> vertices;
        std::set<boost::shared_ptr<Face> > faces;
        const int ID; 
        Edge() : ID(STATIC_ID++) { };
    };

    class Face    
    {
    private:
        static int STATIC_ID;
    public:
        std::set<boost::shared_ptr<Edge> > edges;
        const int ID;
        Face() : ID(STATIC_ID++) { };
    };

    // D für Parameter-Dimension, K für Grad
    // TODO Geometrie ebenfalls als Template-Parameter
    template <unsigned int D = 2, unsigned int K = 3>
    class SimpleTMesh : public Mesh
    {
    private:
        std::set<boost::shared_ptr<Vertex<> > > m_vertices;
        const unsigned int m_K;
    public:
        std::deque<std::pair<IPoint, PPoint> > m_insertionPoints;        // Friend?

    public:

        // Standard-Konstruktor
        SimpleTMesh(void);
        SimpleTMesh(std::set<boost::shared_ptr<Vertex<> > > vertices);
        SimpleTMesh(const Grid& grid);

        // Copy-Konstruktor
        SimpleTMesh(const SimpleTMesh& tmesh);
        // Destruktor
        ~SimpleTMesh(void);
        // Assignment operator
        SimpleTMesh& operator=(const SimpleTMesh& SimpleTMesh);

        /**
        * Erzeugt ein vollbesetzes Kontrollnetz der Größe cols (Spalten) x rows (Zeilen) aus dem Punktegitter. Falls der vector points leer ist, dann werden neue Punkte erzeugt
        * \param cols Spalten
        * \param rows Zeilen
        * \param points Punktegitter
        */
        SimpleTMesh(const std::vector<GPoint>& points, boost::array<unsigned int, D> size);

        /**
        * Setzt eine Kopie des Kontrollpunkts an der Bildkoordinate uv.
        * \param imPt Bildkoordinate
        * \param p Kontrollpunkt
        */
        //void setPoint(Point2d imPt, const GPoint& p);


        /**
        * \brief Lokale Verfeinerung
        * Local Refinement am Bildpunkt imPt.
        * Achtung: Der Bildpunkt markiert lediglich den ungefähren Ort der lokalen Verfeinerung. 
        * Die genaue Refine-Koordinate liegt jeweils mittig auf bzw in dem verfeinerten Objekt.
        * \param iPt Bildbereich-Koordinate
        * \return Liste mit Vertices mit ggf unbestimmten Geometrie-Koordinaten
        */
        void refineByImage(const IPoint& iPt);

        /**
        * \brief Lokale Verfeinerung
        * Local Refinement am Parameterpunkt pPt.
        * \param pPt Parameterraum-Koordinate
        * \return Liste mit Vertices mit ggf unbestimmten Geometrie-Koordinaten
        */
        //void refineByParameter(Point2d pPt);


        /**
        * \brief Lokale Verfeinerung mit bekannten Bild- und Parameterkoordinaten
        * Local Refinement am Parameterpunkt pPt.
        * \param iPt Bildbereich-Koordinate
        * \param pPt Parameterraum-Koordinate
        */
        //void refine(void);
        //void refine(Point2d pPt,Point2d iPt);

        /**
        * \brief  Datei-Import
        * Erzeugt ein T-Mesh aus einer Textdatei
        * \param filename Dateiname
        * \return T-Mesh
        */
        static boost::shared_ptr<SimpleTMesh> fromFile(const char* filename);

        /**
        * \brief  Datei-Export
        * Speichert das T-Grid als Textdatei
        * \param       filename Dateiname
        */
        void toFile(const char* filename) const;

        /**
        * \brief Aktualisiert die Knotenwerte
        * Die lokalen Knotenwerte werden neu aus dem vorliegenden T-Mesh hergeleitet.
        */
        void updateKnotValues(void);

        /**
        * Aktueller Grad des SimpleTMesh
        * \return K
        */
        unsigned int getK() const;

        void clear();

        std::vector<Polygon> toPolygons(void) const;

        MinMaxG minmaxG(void) const;

        void transformG(const AffTransformG& T);

        MinMaxI minmaxI(void) const;

        void transformI(const AffTransformI& T);

        MinMaxP minmaxP(void) const;

        void transformP(const AffTransformP& T);

        TKnotsN getKnotValues(const boost::shared_ptr<Vertex<> > v) const;

        /**
        * \brief Prüft die Konsistenz der Knotenwerte
        * Testet, ob die Knotenwerte an jedem Vertex mit den Bedingungen des T-Mesh übereinstimmen
        */
        bool checkKnotValueConsistency(void) const;

        /**
        * Methode zum Verbinden gegenüberliegender Vertices
        * (ggf Bildkoordinaten mittig verschieben, oder identisch mit ParamKoord)
        */
        //void insertMissingEdges(const Point2d iPt, const Point2d pPt);
        //void insertMissingEdges(void);

        /**
        * \brief Suchen einer adjazenten Kante
        * Sucht die ausgehende Halbkante an einem Vertex in die gegebene Richtung
        * \return Halbkante
        */
        boost::shared_ptr<Edge> findOutgoingEdge(const boost::shared_ptr<Vertex<> > v, const unsigned int inDirection) const;

        /**
        * \brief Suchen des nächsten Vertex
        * Sucht den folgenden Vertex in die gegebene Richtung
        * \return Vertex
        */
        boost::shared_ptr<Vertex<> > findNextVertex(const boost::shared_ptr<Vertex<> > v, const unsigned int inDirection) const;


        /**
        * \brief Parameterkoordinaten -> Bildkoordinaten
        * Berechnet die Bildkoordinaten (i,j) zu gegebenen Parameterkoordinaten p=(s,t).
        * Diese sind u.U. nicht eindeutig, und daher wird eine Menge an Bildpunkten ermittelt.
        * \return Bildkoordinaten
        */
        IPoint imagePoint(const PPoint& p) const;

        /**
        * \brief Bildkoordinaten -> Parameterkoordinaten
        * Berechnet die Parameterkoordinaten (s,t) zu gegebenen Bildkoordinaten i=(u,v).
        * \return Parameterkoordinaten
        */
        PPoint parameterPoint(const IPoint& i) const;

        // *** neue Methoden ***

        /**
        * Gibt den Kontrollpunkt an der Stelle uv im SimpleTMesh aus.
        * \param imPt Bildkoordinate
        * \return Kontrollpunkt
        */
        boost::shared_ptr<Vertex<> > findVertex(const IPoint& imPt);

        /**
        * Gibt den konstanten Kontrollpunkt an der Stelle uv im SimpleTMesh aus.
        * \param imPt Bildkoordinate
        * \return Kontrollpunkt
        */
        const boost::shared_ptr<Vertex<> > findVertexConst(const IPoint& imPt) const;


        std::set<boost::shared_ptr<Vertex<> > > getVertices();
        std::set<boost::shared_ptr<Edge> > getEdges();
        std::set<boost::shared_ptr<Face> > getFaces();
        std::set<const boost::shared_ptr<Vertex<> > > getVerticesConst() const;
        std::set<const boost::shared_ptr<Edge> > getEdgesConst() const;
        std::set<const boost::shared_ptr<Face> > getFacesConst() const;

        static std::pair<MinMaxI,MinMaxP> minmax(const boost::shared_ptr<Edge>);
        static std::pair<MinMaxI,MinMaxP> minmax(const boost::shared_ptr<Face>);

        static boost::shared_ptr<SimpleTMesh<> > getSampleTMesh();
        static boost::shared_ptr<SimpleTMesh<> > getSampleTMesh2();
    };

}

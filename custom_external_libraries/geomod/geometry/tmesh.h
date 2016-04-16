#pragma once

#include "../geomodcore.h"
#include "mesh.h"
#include "grid.h"
#include "tmesh_definition.h"
#include "tmesh/tmeshedgerefine.h"
#include "tmesh/tmeshfacerefine.h"
#include "tmesh/knotvaluesdelegator.h"
#include "tmesh/tmeshsamples.h"

namespace geomodcore {

	enum Direction {
		NORTH, SOUTH, EAST, WEST
	};

	class TMesh : public Mesh, public TMeshArrangement
	{

	private:
		static int count;
		// Point-Location-Objekt
		NaiveTMeshPointLocation pl;

		// Grad des TMesh (wird für Knotenwert-Herleitung benötigt
		unsigned int K;

	public:

		std::deque<std::pair<PPoint2d, PPoint2d> > insertionPoints;

		// Standard-Konstruktor
		TMesh(void);
		// Copy-Konstruktor
		TMesh(const TMesh& tmesh);
		// Destruktor
		~TMesh(void);
		// Assignment operator
		TMesh& operator=(const TMesh& TMesh);

		/**
		* Erzeugt ein vollbesetzes Kontrollnetz der Größe cols (Spalten) x rows (Zeilen) aus dem Punktegitter. Falls der vector points leer ist, dann werden neue Punkte erzeugt
		* \param cols Spalten
		* \param rows Zeilen
		* \param points Punktegitter
		*/
		TMesh(const std::vector<GPoint*>& points, const unsigned int rows, const unsigned int cols);

		/**
		* Setzt eine Kopie des Kontrollpunkts an der Bildkoordinate uv.
		* \param imPt Bildkoordinate
		* \param p Kontrollpunkt
		*/
		void setPoint(PPoint2d imPt, const GPoint& p);

		/**
		* Gibt den konstanten Kontrollpunkt an der Stelle uv im TMesh aus.
		* \param imPt Bildkoordinate
		* \return Kontrollpunkt
		*/
		const GPoint* const getPointConstRef(PPoint2d imPt) const;

		/**
		* Gibt den Kontrollpunkt an der Stelle uv im TMesh aus.
		* \param imPt Bildkoordinate
		* \return Kontrollpunkt
		*/
		GPoint* getPoint(PPoint2d imPt);

		/**
		* Gibt die Kontrollpunkte des TMesh aus.
		* \return Kontrollpunkte
		*/
		const std::vector<GPoint*>* const getPointsConstRef() const;


		/**
		* \brief Lokale Verfeinerung
		* Local Refinement am Bildpunkt imPt.
		* Achtung: Der Bildpunkt markiert lediglich den ungefähren Ort der lokalen Verfeinerung. 
		* Die genaue Refine-Koordinate liegt jeweils mittig auf bzw in dem verfeinerten Objekt.
		* \param iPt Bildbereich-Koordinate
		* \return Liste mit Vertices mit ggf unbestimmten Geometrie-Koordinaten
		*/
		void refineByImage(PPoint2d iPt);

		/**
		* \brief Lokale Verfeinerung
		* Local Refinement am Parameterpunkt pPt.
		* \param pPt Parameterraum-Koordinate
		* \return Liste mit Vertices mit ggf unbestimmten Geometrie-Koordinaten
		*/
		void refineByParameter(PPoint2d pPt);


		/**
		* \brief Lokale Verfeinerung mit bekannten Bild- und Parameterkoordinaten
		* Local Refinement am Parameterpunkt pPt.
		* \param iPt Bildbereich-Koordinate
		* \param pPt Parameterraum-Koordinate
		*/
		void refine(void);
		void refine(PPoint2d pPt,PPoint2d iPt);

		/**
		* \brief  Datei-Import
		* Erzeugt ein T-Mesh aus einer Textdatei
		* \param filename Dateiname
		* \return T-Mesh
		*/
		static TMesh* fromFile(const char* const filename);

		/**
		* \brief  Datei-Export
		* Speichert das T-Grid als Textdatei
		* \param       filename Dateiname
		*/
		void toFile(const char* const filename) const;

		/**
		* \brief Aktualisiert die Knotenwerte
		* Die lokalen Knotenwerte werden neu aus dem vorliegenden T-Mesh hergeleitet.
		*/
		void updateKnotValues(void);

		/**
		* Aktueller Grad des TMesh
		* \return K
		*/
		unsigned int getK() const;

		/**
		* Ändert Grad des TMesh
		* \param K
		*/
		void setK(unsigned int);

		void clear();

		std::vector<Polygon>* toPolygons(void) const;

		MinMaxXYZ minmaxXYZ(void) const;

		void transformXYZ(const AffTransform3d& T);

		MinMaxI minmaxUV(void) const;

		void transformUV(const AffTransform2d& T);

		MinMaxP minmaxST(void) const;

		void transformST(const AffTransform2d& T);

		static TMesh* getSampleTMesh();
		static TMesh* getSampleTMesh2();
		static TMesh* getSampleTMeshMoscow();

		/**
		* \brief Edge-Refinement
		* Fügt einen neuen Bildpunkt i mit Parameterkoordinaten p auf der übergebenen Kante ein
		* Es muss sichergestellt werden, dass der Parameterpunkt innerhalb des gleichen TMesh-Objekts wie i liegt!
		*/
		void refineEdge(TMeshArrangement::Halfedge* const edge, const PPoint2d iPt, const PPoint2d pPt);

		/**
		* \brief Face-Refinement
		* Fügt einen neuen Bildpunkt i mit Parameterkoordinaten ein
		*/
		void refineFace(TMeshArrangement::Face* const face, const PPoint2d iPt, const PPoint2d pPt);

		/**
		* \brief Herleitung der Knotenvektoren aus dem TMesh
		* Bestimmt die lokalen horizontalen & vertikalen Knotenwerte aus dem aktuellen Zustand des TMesh am Vertex v.
		* \return horizontale & vertikale Knotenwerte
		*/
		KnotsN getKnotValues(const TMeshArrangement::Vertex* const v) const;
		KnotsN getKnotValues(const PPoint2d i) const;

		/**
		* \brief Prüft die Konsistenz der Knotenwerte
		* Testet, ob die Knotenwerte an jedem Vertex mit den Bedingungen des T-Mesh übereinstimmen
		*/
		bool checkKnotValueConsistency(void) const;

		/**
		* Methode zum Verbinden gegenüberliegender Vertices
		* (ggf Bildkoordinaten mittig verschieben, oder identisch mit ParamKoord)
		*/
		void insertMissingEdges(const PPoint2d iPt, const PPoint2d pPt);
		void insertMissingEdges(void);

		/**
		* \brief Suchen einer adjazenten Kante
		* Sucht die ausgehende Halbkante an einem Vertex in die gegebene Richtung
		* \return Halbkante
		*/
		TMeshArrangement::Halfedge* findEdge(const TMeshArrangement::Vertex* const v, const Direction direction) const;

		/**
		* \brief Suchen des nächsten Vertex
		* Sucht den folgenden Vertex in die gegebene Richtung
		* \return Vertex
		*/
		TMeshArrangement::Vertex* findVertex(const TMeshArrangement::Vertex* const v, const Direction direction) const;

		/**
		* \brief Suche k Vertices an alle Richtungen
		* Sucht die Menge der k benachbarten Vertices in allen Richtungen (inkl. v), sofern diese existieren.
		* \return benachbarte Vertices
		*/
		std::set<TMeshArrangement::Vertex*> findNeighbourVertices(const TMeshArrangement::Vertex* const v, const unsigned int k) const;

		/**
		* \brief Parameterkoordinaten -> Bildkoordinaten
		* Berechnet die Bildkoordinaten (i,j) zu gegebenen Parameterkoordinaten p=(s,t).
		* Diese sind u.U. nicht eindeutig, und daher wird eine Menge an Bildpunkten ermittelt.
		* \return Bildkoordinaten
		*/
		std::set<PPoint2d> imagePoint(const PPoint2d& p) const;

		/**
		* \brief Bildkoordinaten -> Parameterkoordinaten
		* Berechnet die Parameterkoordinaten (s,t) zu gegebenen Bildkoordinaten i=(u,v).
		* \return Parameterkoordinaten
		*/
		PPoint2d parameterPoint(const PPoint2d& i) const;
		/**
		* \brief Minimale/maximale Bild- und Parameterkoordinaten
		* Gibt die minimalen und maximalen Bild- und Parameterkoordinaten der Zelle aus
		* \return Paar der minimalen und maximalen Koordinaten
		*/
		static std::pair<MinMaxI,MinMaxP> minmax(const TMeshArrangement::Halfedge* const);
		static std::pair<MinMaxI,MinMaxP> minmax(const TMeshArrangement::Face* const);

		/**
		* Skaliert die geometrischGeometrie des TMesh auf den Bereich (1,1,1)
		*/
		void scaleToSimplex();

	private:
		/**
		* Weist jedem TMesh-Objekt eine neue ID zu.
		*/
		void reassignIds(void);

	};

}

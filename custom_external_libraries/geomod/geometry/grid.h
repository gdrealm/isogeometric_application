#pragma once

#include "../geomodcore.h"
#include "mesh.h"

namespace geomodcore {

    class Grid : public Mesh {

    private:
        unsigned int cols;
        unsigned int rows;

    protected:

        /**
        * Punkte im Gitter in row-major-Ordnung. Da das Gitter u.U. Lücken enthält, sind die Kontrollpunkte optional (vgl Boost optional)
        */
        std::vector<boost::shared_ptr<GPoint> > m_points;

    public:

        // Standard-Konstruktor
        Grid(void);
        // Copy-Konstruktor
        Grid(const Grid& grid);
        // Destruktor
        ~Grid(void);
        // Assignment operator
        Grid& operator=(const Grid& grid);

        /**
        * Erzeugt ein Kontrollnetz der Größe cols (Spalten) x rows (Zeilen).
        * \param cols Spalten
        * \param rows Zeilen
        */
        Grid(const unsigned int cols, const unsigned int rows);


        /**
        * Erzeugt ein Kontrollnetz aus dem Punktegitter.
        * \param cols Spalten
        * \param rows Zeilen
        * \param points Punktegitter
        */
        Grid(const std::vector<boost::shared_ptr<GPoint> >& points, const unsigned int rows, const unsigned int cols);

        /**
        * Setzt eine Kopie des Kontrollpunkts an der Stelle (i,j).
        * \param col Spalte
        * \param row Zeile
        * \param p Kontrollpunkt
        */
        void setPoint(const unsigned int col, const unsigned int row, const GPoint& p);

        /**
        * Gibt den Kontrollpunkt an der Stelle (i,j) im Raster aus.
        * \param col Spalte
        * \param row Zeile
        * \return Kontrollpunkt
        */
        const boost::shared_ptr<GPoint> getControlPointConstRef(const unsigned int col, const unsigned int row) const;

        /**
        * Gibt den i-ten Kontrollpunkt im Raster aus.
        * \param i Index in Reihen-Spalten-Sortierung
        * \return Kontrollpunkt
        */
        const boost::shared_ptr<GPoint> getControlPointConstRef(const unsigned int i) const;

        
        /**
        * Gibt den Kontrollpunkt an der Stelle (i,j) im Raster aus.
        * \param col Spalte
        * \param row Zeile
        * \return Kontrollpunkt
        */
        boost::shared_ptr<GPoint> getControlPoint(const unsigned int col, const unsigned int row);

        /**
        * Gibt den i-ten Kontrollpunkt im Raster aus.
        * \param i Index in Reihen-Spalten-Sortierung
        * \return Kontrollpunkt
        */
        boost::shared_ptr<GPoint> getControlPoint(const unsigned int i);


        /**
        * Gibt die Kontrollpunkte des Rasters aus.
        * \return Kontrollpunkte
        */
        const std::vector<boost::shared_ptr<GPoint> >& getControlPointsConstRef() const;

        /**
        * Gibt das Kontrollpolygon der j-ten Reihe aus.
        * \param row Reihenindex
        * \return Kontrollpolygon
        */
        Polygon getRowPolygon(const unsigned int row) const;

        /**
        * Gibt das Kontrollpolygon der i-ten Spalte aus.
        * \param col Spaltenindex
        * \return Kontrollpolygon
        */
        Polygon getColPolygon(const unsigned int col) const;

        /**
        * Anzahl der Spalten
        * \return Anzahl der Spalten
        */
        unsigned int getCols() const;

        /**
        * Anzahl der Reihen
        * \return Anzahl der Reihen
        */
        unsigned int getRows() const;

        /**
        * \brief ESRI-ASCII-Grid-Import
        * Importiert ein ESRI-ASCII-Grid als Grid
        * \param filename Dateiname
        * \param noValueData Höhenwert bei NOVALUE_DATA
        */
        static boost::shared_ptr<Grid> readGridFromESRIFile(const std::string filename, const double noValueData = 0);

        std::vector<Polygon> toPolygons(void) const;

        void clear(void);

        MinMaxG minmaxG(void) const;

        void transformG(const AffTransformG& T);

        MinMaxI minmaxI(void) const;

        void transformI(const AffTransformP& T);

        static boost::shared_ptr<Grid> getSampleGrid();
    };

}

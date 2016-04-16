#include "grid.h"

using namespace std;
using namespace geomodcore;

Grid::Grid(void)
{
}

Grid::Grid(const unsigned int cols, const unsigned int rows) 
	: cols(cols), rows(rows)
{
}


Grid::Grid(const std::vector<boost::shared_ptr<GPoint> >& points, const unsigned int rows, const unsigned int cols) 
	: cols(cols), rows(rows)
{
	this->m_points.assign(points.begin(), points.end());
}

Grid::Grid(const Grid& grid)
	: cols(grid.cols), rows(grid.rows)
{
	this->m_points.assign(grid.m_points.begin(), grid.m_points.end());
}

Grid& Grid::operator=(const Grid& grid) 
{
    this->m_points.assign(grid.m_points.begin(), grid.m_points.end());
	return *this;
}

Grid::~Grid(void) {
	clear();
}

unsigned int Grid::getRows() const 
{
	return rows;
}
unsigned int Grid::getCols() const
{
	return cols;
}

void Grid::setPoint(const unsigned int col, const unsigned int row, const GPoint& p) 
{
	m_points[row*cols + col] = boost::shared_ptr<GPoint>(new GPoint(p));
}

boost::shared_ptr<GPoint> Grid::getControlPoint(const unsigned int i)  
{
	return m_points[i];
}

boost::shared_ptr<GPoint> Grid::getControlPoint(const unsigned int col, const unsigned int row) 
{
	return m_points[row*cols+ col];
}


const boost::shared_ptr<GPoint> Grid::getControlPointConstRef(const unsigned int i) const 
{
	return m_points[i];
}

const boost::shared_ptr<GPoint> Grid::getControlPointConstRef(const unsigned int col, const unsigned int row) const 
{
	return m_points[row*cols+ col];
}

const std::vector<boost::shared_ptr<GPoint> >& Grid::getControlPointsConstRef() const
{
	return m_points;
}

Polygon Grid::getRowPolygon(const unsigned int row) const
{
	Polygon rowPoly;
	for (unsigned int col = 0; col < cols; col++) {
		if (!(m_points[row*cols + col]))
		{
			rowPoly.push_back(*m_points[row*cols + col]);
		}
	}
	return Polygon(rowPoly);
}

Polygon Grid::getColPolygon(const unsigned int col) const
{
	Polygon colPoly;
	for (unsigned int row = 0; row < rows; row++) {
		if (!(m_points[row*cols + col]))
		{
			colPoly.push_back(*m_points[row*cols + col]);
		}
	}
	return Polygon(colPoly);
}

std::vector<Polygon> Grid::toPolygons(void) const
{
	assert((cols > 1) && (rows > 1));
	std::vector<Polygon> polygons;

	for (unsigned int row = 0; row < rows - 1; row++) {
		for (unsigned int col = 0; col < cols - 1; col++) {
			Polygon p;
			if (!(getControlPointConstRef(col, row))) p.push_back(*getControlPointConstRef(col, row));
			if (!(getControlPointConstRef(col+1, row))) p.push_back(*getControlPointConstRef(col+1, row));
			if (!(getControlPointConstRef(col+1, row+1))) p.push_back(*getControlPointConstRef(col+1, row+1));
			if (!(getControlPointConstRef(col, row+1))) p.push_back(*getControlPointConstRef(col, row+1));
			polygons.push_back(p);
		}
	}

	return polygons;
}

void Grid::clear(void) 
{
	m_points.clear();
	rows = 0;
	cols = 0;
}

MinMaxG Grid::minmaxG() const 
{
	GPoint minG = {_MAXDOUBLE_, _MAXDOUBLE_, _MAXDOUBLE_};
	GPoint maxG = {_MINDOUBLE_, _MINDOUBLE_, _MINDOUBLE_};
	for (unsigned int i = 0; i < m_points.size(); ++i) {
		if (!(m_points[i])) {
			const GPoint& p = *m_points[i];
			minG[0] = min(minG[0], p[0]);
			minG[1] = min(minG[1], p[1]);
			minG[2] = min(minG[2], p[2]);
			maxG[0] = max(maxG[0], p[0]);
			maxG[1] = max(maxG[1], p[1]);
			maxG[2] = max(maxG[2], p[2]);
		}
	}
	return MinMaxG(minG,maxG);
}

void Grid::transformG(const AffTransformG& T) 
{
	for (unsigned int i = 0; i < m_points.size(); ++i) {
		if (!(m_points[i])) {
		    GPoint& p = *m_points[i];
		    p = transform(p, T);
		}
	}
}

MinMaxI Grid::minmaxI(void) const 
{
	PPoint min = {0,0};
	PPoint max = {this->cols, this->rows};
	return MinMaxI(min, max);
}

void Grid::transformI(const AffTransformP& T) 
{
	cout << "2D-Transformation auf Grid-Netze werden nicht unterstützt!" << endl;
}

boost::shared_ptr<Grid> Grid::getSampleGrid() 
{
	const int n = GRID_SAMPLE_SIZE;
	std::vector<boost::shared_ptr<GPoint> > points;
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
            boost::shared_ptr<GPoint> p = boost::shared_ptr<GPoint>(new GPoint());
            (*p)[0] = col;
            (*p)[1] = col;
            (*p)[2] = std::cos((double)col) * std::sin((double)row);
			points.push_back(p);
		}
	}
	return boost::shared_ptr<Grid>(new Grid(points, n, n));
}

boost::shared_ptr<Grid> Grid::readGridFromESRIFile(const std::string filename, const double noValueData)
{
	cout << "Importiere SimpleTMesh aus " << filename << endl;

	// Arbeitsvariablen
	std::fstream file;
	file.open(filename.c_str());
	string line;
	std::vector<std::string> linesplit;

	// Datenvariablen
	unsigned int cols = 0;
	unsigned int rows = 0;
	double x0 = 0;
	double y0 = 0;
	double cellsize = 0;
	std::string noData = "";
	std::vector<boost::shared_ptr<GPoint> > points;

	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" ,\t"));
	cols = boost::lexical_cast<int>(linesplit[1]);
	cout << "> Spalten: " << cols << endl;

	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	rows = boost::lexical_cast<int>(linesplit[1]);
	cout << "> Reihen: " << rows << endl;
	cout << "ÄNDERE ROWS AB!" << endl;
	rows = 20;
	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	if (linesplit[1].find(",") < linesplit[1].size())
	{
		linesplit[1].replace(linesplit[1].find(","),1,".");
	}
	x0 = boost::lexical_cast<double>(linesplit[1]);
	cout << "> Minimale x-Koordinate: " << x0 << endl;

	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	if (linesplit[1].find(",") < linesplit[1].size())
	{
		linesplit[1].replace(linesplit[1].find(","),1,".");
	}
	y0 = boost::lexical_cast<double>(linesplit[1]);
	cout << "> Minimale y-Koordinate: " << y0 << endl;

	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	if (linesplit[1].find(",") < linesplit[1].size())
	{
		linesplit[1].replace(linesplit[1].find(","),1,".");
	}
	cellsize = boost::lexical_cast<double>(linesplit[1]);
	cout << "> Cell-Size: " << cellsize << endl;

	do 
	{
		getline(file, line);
	} while (line.empty());
	boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
	noData = linesplit[1];
	cout << "> NODATA_value: " << noData << endl;

	// Vorbereitung
	points.assign(cols*rows, boost::shared_ptr<GPoint>(new GPoint()));

	vector<string> lines;
	for (unsigned int currentRow = 0; currentRow < rows; currentRow++)
	{
		getline(file, line);
		lines.push_back(line);
	}

	cout << "Importiere ";		
	//#pragma omp parallel for
	for (unsigned int currentRow = 0; currentRow < rows; currentRow++)
	{
		const string& line = lines[currentRow];
		boost::algorithm::split(linesplit, line, boost::is_any_of(" "));
		for (unsigned int currentCol = 0; currentCol < cols; currentCol++)
		{
			double val;
			if (linesplit[currentCol] == noData)
			{
				val = noValueData;
			} else {
				if (linesplit[currentCol].find(",") < linesplit[currentCol].size())
				{
					linesplit[currentCol].replace(linesplit[currentCol].find(","),1,".");
				}
				val = boost::lexical_cast<double>(linesplit[currentCol]);
			}
			boost::shared_ptr<GPoint> p = boost::shared_ptr<GPoint>(new GPoint());
            (*p)[0] = x0 + currentCol * cellsize;
            (*p)[1] = y0 + currentRow * cellsize;
            (*p)[2] = val;
			points[currentRow * cols + currentCol].swap(p);
		}
		cout << ".";
	}
	cout << " fertig" << endl;

	assert(points.size() == rows*cols);
	return boost::shared_ptr<Grid>(new Grid(points, rows, cols));
}

#include "tmeshsamples.h"

using namespace std;
using namespace geomodcore;


SimpleTMesh<>* TMeshSamples::getSampleTMesh()
{
	const Grid* const sampleGrid = Grid::getSampleGrid();
	return new SimpleTMesh<>(sampleGrid);
};


SimpleTMesh<>* TMeshSamples::getSampleTMesh2()
{
	string filename;
	filename.append("C:/Users/").append(getenv("USERNAME")).append("/Projekte/TSplines/tmesh_sample_paper.txt");
	return SimpleTMesh<>::fromFile(filename.data());
};

/*

SimpleTMesh<>* TMeshSamples::getSampleTMeshMoscow()
{
	const unsigned int n = 100;
	std::vector<GPoint*> points;
	points.assign(n*n, new GPoint());


	for (int row = 0; row < 50; row++) {
		if (row < 25) {
			for (int col = 15; col < 35; col++) {
				points[col*n+row] = new GPoint(col + (row - 25) / 3., row, -5 + abs(25 - col) * 0.5);
			}
		} else {
			for (int col = 0; col < 15; col++) {
				if (row < 30) {
					points[col*n+row] = new GPoint(col, row, 25 - row);
				} else if (row < 45) {
					points[col*n+row] = new GPoint(col, row, -5.);
				} else {
					points[col*n+row] = new GPoint(col, row, row - 50);
				}
			}
			for (int col = 35; col < 50; col++) {
				if (row < 30) {
					points[col*n+row] = new GPoint(col, row, 25 - row);
				} else if (row < 45) {
					points[col*n+row] = new GPoint(col, row, -5.);
				} else {
					points[col*n+row] = new GPoint(col, row, row - 50);
				}
			}
			for (int col = 15; col < 35; col++) {
				if (row < 30) {
					if (col < 25) {
						double d0 = - (col-15) / 2.;
						double e0 = d0 * (30 - row) / 5. - 5. * (row - 25) / 5.;
						points[col*n+row] = new GPoint(col, row, e0);
					} else if (col < 35) {
						double d0 = - (35-col) / 2.;
						double e0 = d0 * (30 - row) / 5. - 5. * (row - 25) / 5.;
						points[col*n+row] = new GPoint(col, row, e0);
					}
				} else if (row < 45) {
					points[col*n+row] = new GPoint(col, row, -5.);
				} else {
					points[col*n+row] = new GPoint(col, row, row - 50);
				}
			}
		}
	}
	for (unsigned int i = 0; i < n*n; i++) {
		if (points[i] != nullptr) {
			points[i] = new GPoint(points[i]->x(), points[i]->y(), points[i]->z() + 0.3 * (std::rand() / RAND_MAX));
		}
	}
	return new SimpleTMesh<>(points, n, n);
}

*/
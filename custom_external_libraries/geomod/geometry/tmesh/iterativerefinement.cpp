#include "iterativerefinement.h"

using namespace std;
using namespace geomodcore;

IterativeRefinement::IterativeRefinement(const unsigned int maxRefineDepth) : maxDepth(maxRefineDepth)
{
}

SimpleTMesh<>* IterativeRefinement::refine(const BSplineSurface& bspline)
{
	cout << "Diskretisierung der B-Spline-Fläche ... ";
	const Grid* const bsplineSurfaceDisc = bspline.getDiscretization(REFINEMENT_DISCRETIZATION_SIZE, REFINEMENT_DISCRETIZATION_SIZE);
	cout << " fertig." << endl;
	const unsigned int cols = bspline.getGridConstRef().getCols();
	const unsigned int rows = bspline.getGridConstRef().getRows();

	vector<GPoint*> points;
	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(0,0)));
	//	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef((cols-1)/2,0)));
	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(cols-1,0)));

	//	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(0,(rows-1)/2)));
	//	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef((cols-1)/2,(rows-1)/2)));
	//	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(cols-1,(rows-1)/2)));

	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(0,rows-1)));
	//	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef((cols-1)/2,rows-1)));
	points.push_back(new GPoint(*bspline.getGridConstRef().getControlPointConstRef(cols-1,rows-1)));

	SimpleTMesh<>* tmesh = new SimpleTMesh<>(bspline.getGridConstRef());
	const AffTransformP T = { cols, 0, 0, 0, rows, 0 };
	tmesh->transformI(T);
	// TODO
	cout << "Speichere SimpleTMesh nach Transformation" << endl;
	tmesh->toFile("C:/Users/asche/Desktop/refine.txt");
	TSplineSurface<> tspline(*tmesh);

	double maxError = MAXDOUBLE;
	unsigned int* x = new unsigned int;
	unsigned int* y = new unsigned int;

	for (unsigned int d = 0; d <= maxDepth; d++)
	{
		const unsigned int n = 1 << d;
		cout << "n = " << n << endl;
		for (unsigned int i = 0; i < n * n; i++)
		{
			HilbertCurve::d2xy(n, i, x, y);
			const unsigned int u = (cols / n * *x) + (cols >> (d+1));
			const unsigned int v = (rows / n * *y) + (rows >> (d+1));
			SimpleTMesh<>* tmeshBeforeRefinement = new SimpleTMesh<>(*tmesh);
			IPoint refineIPt = {u,v};
			tmesh->refineByImage(refineIPt);
			TSplineSurface<> tSplineSurfaceRefined(*tmesh);
			const Grid* const tsplineSurfaceDisc = tSplineSurfaceRefined.getDiscretization(REFINEMENT_DISCRETIZATION_SIZE, REFINEMENT_DISCRETIZATION_SIZE); 
			double newError = diff(*bsplineSurfaceDisc, *tsplineSurfaceDisc);

			if ( newError > maxError ) {
				delete tmesh;
				tmesh = tmeshBeforeRefinement;
			} else {
				cout << "Optimierung: " << newError << endl;
				maxError = newError;
				delete tmeshBeforeRefinement;
				tmesh->toFile("C:/Users/asche/Desktop/refine.txt");
			}
			delete tsplineSurfaceDisc;
			cout << "***" << d << " " << i << endl;
		}
		tmesh->toFile("C:/Users/asche/Desktop/refine.txt");
	}
	delete bsplineSurfaceDisc;
	cout << "Globaler Fehler: " << maxError << endl;

	return tmesh;
}

double IterativeRefinement::diff(const Grid& grid1, const Grid& grid2) const
{
	assert((grid1.getCols() == grid2.getCols()) && (grid1.getRows() == grid2.getRows()));
	const int cols = grid1.getCols();
	const int rows = grid1.getRows();
	double maxError = MINDOUBLE;
	double globalError = 0.0;
	const double xStart = min(grid1.getControlPointConstRef(0, 0)->at(0), grid1.getControlPointConstRef(cols - 1, rows - 1)->at(0));
	const double yStart = min(grid1.getControlPointConstRef(0, 0)->at(1), grid1.getControlPointConstRef(cols - 1, rows - 1)->at(1));
	const double xEnd = max(grid1.getControlPointConstRef(0, 0)->at(0), grid1.getControlPointConstRef(cols - 1, rows - 1)->at(0));
	const double yEnd = max(grid1.getControlPointConstRef(0, 0)->at(1), grid1.getControlPointConstRef(cols - 1, rows - 1)->at(1));
	const double dx = (xEnd - xStart) / cols;
	const double dy = (yEnd - yStart) / rows;
	int count = 0;

	bool hasNullPointsGrid1 = false;
	bool hasNullPointsGrid2 = false;
	for (vector<GPoint*>::const_iterator pIt = grid1.getControlPointsConstRef().cbegin(); pIt != grid1.getControlPointsConstRef().cend(); pIt++)
	{
		hasNullPointsGrid1 = hasNullPointsGrid1 | (*pIt == nullptr);
	}
	for (vector<GPoint*>::const_iterator pIt = grid2.getControlPointsConstRef().cbegin(); pIt != grid2.getControlPointsConstRef().cend(); pIt++)
	{
		hasNullPointsGrid2 = hasNullPointsGrid2 | (*pIt == nullptr);
	}
	if (hasNullPointsGrid1 == false && hasNullPointsGrid2 == false)
	{
		cout << "Berechne schnelle Shephard-Interpolation...";
#pragma omp parallel for reduction(+:globalError), reduction(+:count)
		for (int i = 0; i < (int) grid1.getControlPointsConstRef().size(); i++)
		{
			const GPoint& p = *grid1.getControlPointsConstRef()[i];
			double x = p[0];			
			double y = p[1];
			const double z1 = shepardInterpolationWithoutNulls(grid1, x, y);
			const double z2 = shepardInterpolationWithoutNulls(grid2, x, y);
			const double e = (z1 - z2);
			maxError = max(e * e, maxError);
			globalError += (e * e);
#pragma omp atomic 
			count++;
		}
	} else {
		cout << "Berechne Shephard-Interpolation...";
		for (int i = 0; i < (int) grid1.getControlPointsConstRef().size(); i++)
		{
			const GPoint* const p1 = grid1.getControlPointsConstRef()[i];
			const GPoint* const p2 = grid2.getControlPointsConstRef()[i];

			const double z1 = (p1 != nullptr ? (*p1)[2] : 0);
			const double z2 = (p2 != nullptr ? (*p2)[2] : 0);
			const double e = (z1 - z2);
			maxError = max(e * e, maxError);
			globalError += (e * e);
			count++;
		}
	}

	cout << "fertig." << endl;
	globalError = sqrt(globalError) / count;
	return globalError;
}

double IterativeRefinement::shepardInterpolation(const Grid& grid, double x, double y) const
{
	double z = 0.0;
	double D = 0.0;
	for (vector<GPoint*>::const_iterator pIt = grid.getControlPointsConstRef().cbegin(); pIt != grid.getControlPointsConstRef().cend(); pIt++)
	{
		const GPoint* const p = *pIt;
		if (p != nullptr) {
			const double d = 1.0 / (pow((*p)[0] - x, 4.0) + pow((*p)[1] - y, 4.));
			D += d;
			z += (*p)[2] * d;
		}
	}
	return z / D;
}


double IterativeRefinement::shepardInterpolationWithoutNulls(const Grid& grid, double x, double y) const
{
	double z = 0.0;
	double D = 0.0;
	for (vector<GPoint*>::const_iterator pIt = grid.getControlPointsConstRef().cbegin(); pIt != grid.getControlPointsConstRef().cend(); pIt++)
	{
		const GPoint* const p = *pIt;
		const double d = 1.0 / (pow((*p)[0] - x, 4.0) + pow((*p)[1] - y, 4.));
		D += d;
		z += (*p)[2] * d;
	}
	return z / D;
}


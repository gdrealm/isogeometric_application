#pragma once

#include "../../geomodcore.h"
#include "../../surface/bsplinesurface.h"
#include "../../surface/tsplinesurface.h"
#include "simpletmesh.h"
#include "../tmesh/hilbertcurve.h"

namespace geomodcore {

	class IterativeRefinement
	{
	private:
		const unsigned int maxDepth;

	public:
		IterativeRefinement(const unsigned int maxRefineDepth = 5);

		SimpleTMesh<>* refine(const BSplineSurface& bspline);

	private:
		double diff(const geomodcore::Grid& grid1, const geomodcore::Grid& grid2) const;
		double shepardInterpolation(const geomodcore::Grid& grid, double x, double y) const;
		double shepardInterpolationWithoutNulls(const Grid& grid, double x, double y) const;

	};

}
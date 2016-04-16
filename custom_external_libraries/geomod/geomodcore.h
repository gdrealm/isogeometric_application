#pragma once
#define _USE_MATH_DEFINES

#include "geomodcore_global.h"

#include <cstdio>
#include <cmath>
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <limits>

// #include <array>
#include <boost/array.hpp>
#include <vector>
#include <map>
#include <set>
// #include <unordered_map>
#include <boost/unordered_map.hpp>
#include <deque>
#include <queue>
#include <iterator>

#include "geometry/geometry_definitions.h"

#include <omp.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/tuple/tuple.hpp>

namespace geomodcore {
	
	const double _EPSILON_ = std::numeric_limits<double>::epsilon();
	const double _NAN_ = std::numeric_limits<double>::quiet_NaN();
	const double _INFINITY_ = std::numeric_limits<double>::infinity();
	const double _MINDOUBLE_ = std::numeric_limits<double>::min();
	const double _MAXDOUBLE_ = std::numeric_limits<double>::max();

	const unsigned int DISCRETIZATION_SIZE = 100;
	const unsigned int REFINEMENT_DISCRETIZATION_SIZE = 50;
	const unsigned int GRID_SAMPLE_SIZE = 5;

	class GEOMODCORE_EXPORT GeomodCore
	{
	public:
		GeomodCore();
		~GeomodCore();

	private:

	};
}


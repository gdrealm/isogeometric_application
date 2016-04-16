#pragma once

#include "../../geomodcore.h"

namespace geomodcore {

	class HilbertCurve
	{
	public:
		//convert (x,y) to d
		static unsigned int xy2d (const unsigned int n, const unsigned int x, const unsigned int y);
		//convert d to (x,y)
		static void d2xy(const unsigned int n, const unsigned int d, unsigned int *x, unsigned int *y);

	private:
		//rotate/flip a quadrant appropriately
		static void rot(const unsigned int n, unsigned int *x, unsigned int *y, const unsigned int rx, const unsigned int ry);
	};
}

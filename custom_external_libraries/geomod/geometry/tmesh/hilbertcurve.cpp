#pragma once

#include "hilbertcurve.h"

namespace geomodcore {

	unsigned int HilbertCurve::xy2d (const unsigned int n, unsigned int x, unsigned int y)
	{
		unsigned int rx, ry, s, d=0;
		for (s=n/2; s>0; s/=2) {
			rx = (x & s) > 0;
			ry = (y & s) > 0;
			d += s * s * ((3 * rx) ^ ry);
			rot(s, &x, &y, rx, ry);
		}
		return d;
	}

	void HilbertCurve::d2xy(const unsigned int n, const unsigned int d, unsigned int *x, unsigned int *y)
	{
		unsigned int rx, ry, s, t=d;
		*x = *y = 0;
		for (s=1; s<n; s*=2) {
			rx = 1 & (t/2);
			ry = 1 & (t ^ rx);
			rot(s, x, y, rx, ry);
			*x += s * rx;
			*y += s * ry;
			t /= 4;
		}
	}

	void HilbertCurve::rot(const unsigned int n, unsigned int *x, unsigned int *y, const unsigned int rx, const unsigned int ry)
	{
		if (ry == 0) {
			if (rx == 1) {
				*x = n-1 - *x;
				*y = n-1 - *y;
			}

			//Swap x and y
			int t  = *x;
			*x = *y;
			*y = t;
		}
	}

}

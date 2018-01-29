
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspderiv(d, c, k):

    # Local Variables: tmp, c, nk, d, mc, i, k, nc, dc, dk
    # Function calls: zeros, bspderiv, numel, size
    #% BSPDERIV:  B-Spline derivative.
    #% 
    #%  MATLAB SYNTAX:
    #% 
    #%         [dc,dk] = bspderiv(d,c,k)
    #%  
    #%  INPUT:
    #% 
    #%    d - degree of the B-Spline
    #%    c - control points          double  matrix(mc,nc)
    #%    k - knot sequence           double  vector(nk)
    #% 
    #%  OUTPUT:
    #% 
    #%    dc - control points of the derivative     double  matrix(mc,nc)
    #%    dk - knot sequence of the derivative      double  vector(nk)
    #% 
    #%  Modified version of Algorithm A3.3 from 'The NURBS BOOK' pg98.
    #%
    #%    Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton
    #%
    #%    This program is free software: you can redistribute it and/or modify
    #%    it under the terms of the GNU General Public License as published by
    #%    the Free Software Foundation, either version 3 of the License, or
    #%    (at your option) any later version.
    #%    This program is distributed in the hope that it will be useful,
    #%    but WITHOUT ANY WARRANTY; without even the implied warranty of
    #%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #%    GNU General Public License for more details.
    #%
    #%    You should have received a copy of the GNU General Public License
    #%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    [mc, nc] = matcompat.size(c)
    nk = numel(k)
    #%
    #% int bspderiv(int d, double *c, int mc, int nc, double *k, int nk, double *dc,
    #%              double *dk)
    #% {
    #%   int ierr = 0;
    #%   int i, j, tmp;
    #%
    #%   // control points
    #%   double **ctrl = vec2mat(c,mc,nc);
    #%
    #%   // control points of the derivative
    dc = np.zeros(mc, (nc-1.))
    #%   double **dctrl = vec2mat(dc,mc,nc-1);
    #%
    for i in np.arange(0., (nc-2.)+1):
        #%   for (i = 0; i < nc-1; i++) {
        
    #%     }
    #%   }
    #%
    dk = np.zeros(1., (nk-2.))
    #%   j = 0;
    dk[0:nk-2.] = k[1:nk-1.]
    #%   for (i = 1; i < nk-1; i++)
    #%     dk[j++] = k[i];
    #%
    #%   freevec2mat(dctrl);
    #%   freevec2mat(ctrl);
    #%
    #%   return ierr;
    #% }
    return [dc, dk]
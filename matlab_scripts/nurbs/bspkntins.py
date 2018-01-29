
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspkntins(d, c, k, u):

    # Local Variables: alfa, ii, ik, ind, ic, tmp, nk, nc, nu, C, jj, a, c, b, d, mc, ss, k, m, l, n, r, u
    # Function calls: sort, bspkntins, abs, zeros, numel, findspan, size
    #% BSPKNTINS:  Insert knots into a B-Spline
    #%
    #% Calling Sequence:
    #% 
    #%   [ic,ik] = bspkntins(d,c,k,u)
    #%
    #%  INPUT:
    #% 
    #%    d - spline degree             integer
    #%    c - control points            double  matrix(mc,nc)      
    #%    k - knot sequence             double  vector(nk) 
    #%    u - new knots                 double  vector(nu)               
    #% 
    #%  OUTPUT:
    #% 
    #%    ic - new control points double  matrix(mc,nc+nu) 
    #%    ik - new knot sequence  double  vector(nk+nu)
    #% 
    #%  Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164.
    #% 
    #%    Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton, 2010-2016 Rafael Vazquez
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
    u = np.sort(u)
    nu = numel(u)
    nk = numel(k)
    #% 
    #% int bspkntins(int d, double *c, int mc, int nc, double *k, int nk,
    #%               double *u, int nu, double *ic, double *ik)
    #% {
    #%   int ierr = 0;
    #%   int a, b, r, l, i, j, m, n, s, q, ind;
    #%   double alfa;
    #%
    #%   double **ctrl  = vec2mat(c, mc, nc);
    ic = np.zeros(mc, (nc+nu))
    #%   double **ictrl = vec2mat(ic, mc, nc+nu);
    ik = np.zeros(1., (nk+nu))
    #%
    n = nc-1.
    #%   n = nc - 1;
    r = nu-1.
    #%   r = nu - 1;
    #%
    m = n+d+1.
    #%   m = n + d + 1;
    a = findspan(n, d, u[0], k)
    #%   a = findspan(n, d, u[0], k);
    b = findspan(n, d, u[int((r+1.))-1], k)
    #%   b = findspan(n, d, u[r], k);
    b = b+1.
    #%   ++b;
    #%
    #%   for (q = 0; q < mc; q++)  {
    ic[:,0:a-d+1.] = c[:,0:a-d+1.]
    #%     for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q];
    ic[:,int(b+nu)-1:nc+nu] = c[:,int(b)-1:nc]
    #%     for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q];
    #%   }
    ik[0:a+1.] = k[0:a+1.]
    #%   for (j = 0; j <= a; j++)   ik[j] = k[j];
    ik[int(b+d+nu+1.)-1:m+nu+1.] = k[int(b+d+1.)-1:m+1.]
    #%   for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];
    #%
    ii = b+d-1.
    #%   i = b + d - 1;
    ss = ii+nu
    #%   s = b + d + r;
    for jj in np.arange(r, (0.)+(-1.), -1.):
        #%   for (j = r; j >= 0; j--) {
        
    #%   }
    #%
    #%   freevec2mat(ctrl);
    #%   freevec2mat(ictrl);
    #%
    #%   return ierr;
    #% }
    return [ic, ik, C]
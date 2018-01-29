
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspdegelev(d, c, k, t):

    # Local Variables: alf, inv, lbz, ii, ik, bpts, alfs, ph2, ic, tmp1, kind, tmp2, oldr, cind, nc, tr, mul, ph, save, rbz, numer, bezalfs, Nextbpts, mpi, gam, ebpts, a, c, b, last, d, mc, i, kj, k, j, m, mh, n, q, s, r, t, den, first, ua, bet, ub
    # Function calls: bincoeff, min, max, floor, zeros, bspdegelev, size
    #% BSPDEGELEV:  Degree elevate a univariate B-Spline. 
    #% 
    #% Calling Sequence:
    #% 
    #%   [ic,ik] = bspdegelev(d,c,k,t)
    #% 
    #%   INPUT:
    #% 
    #%   d - Degree of the B-Spline.
    #%   c - Control points, matrix of size (dim,nc).
    #%   k - Knot sequence, row vector of size nk.
    #%   t - Raise the B-Spline degree t times.
    #% 
    #%   OUTPUT:
    #%
    #%   ic - Control points of the new B-Spline. 
    #%   ik - Knot vector of the new B-Spline.
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
    #%
    #% int bspdegelev(int d, double *c, int mc, int nc, double *k, int nk,
    #%                int t, int *nh, double *ic, double *ik)
    #% {
    #%   int row,col;
    #%
    #%   int ierr = 0;
    #%   int i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, cind, oldr, mul;
    #%   int n, lbz, rbz, save, tr, kj, first, kind, last, bet, ii;
    #%   double inv, ua, ub, numer, den, alf, gam;
    #%   double **bezalfs, **bpts, **ebpts, **Nextbpts, *alfs;
    #%
    #%   double **ctrl  = vec2mat(c, mc, nc);
    #% ic = zeros(mc,nc*(t));                                  %   double **ictrl = vec2mat(ic, mc, nc*(t+1));
    #%
    n = nc-1.
    #%   n = nc - 1;
    #%
    bezalfs = np.zeros((d+1.), (d+t+1.))
    #%   bezalfs = matrix(d+1,d+t+1);
    bpts = np.zeros(mc, (d+1.))
    #%   bpts = matrix(mc,d+1);
    ebpts = np.zeros(mc, (d+t+1.))
    #%   ebpts = matrix(mc,d+t+1);
    Nextbpts = np.zeros(mc, (d+1.))
    #%   Nextbpts = matrix(mc,d+1);
    alfs = np.zeros(d, 1.)
    #%   alfs = (double *) mxMalloc(d*sizeof(double));
    #%
    m = n+d+1.
    #%   m = n + d + 1;
    ph = d+t
    #%   ph = d + t;
    ph2 = np.floor((ph/2.))
    #%   ph2 = ph / 2;
    #%
    #%   // compute bezier degree elevation coefficeients
    bezalfs[0,0] = 1.
    #%   bezalfs[0][0] = bezalfs[ph][d] = 1.0;
    bezalfs[int((d+1.))-1,int((ph+1.))-1] = 1.
    #%
    for i in np.arange(1., (ph2)+1):
        #%   for (i = 1; i <= ph2; i++) {
        
    #%   }
    #%
    for i in np.arange(ph2+1., (ph-1.)+1):
        #%   for (i = ph2+1; i <= ph-1; i++) {
        
    #%   }
    #%
    mh = ph
    #%   mh = ph;      
    kind = ph+1.
    #%   kind = ph+1;
    r = -1.
    #%   r = -1;
    a = d
    #%   a = d;
    b = d+1.
    #%   b = d+1;
    cind = 1.
    #%   cind = 1;
    ua = k[0]
    #%   ua = k[0]; 
    #%
    for ii in np.arange(0., (mc-1.)+1):
        #%   for (ii = 0; ii < mc; ii++)
        
    #%
    for i in np.arange(0., (ph)+1):
        #%   for (i = 0; i <= ph; i++)
        
    #%
    #%   // initialise first bezier seg
    for i in np.arange(0., (d)+1):
        #%   for (i = 0; i <= d; i++)
        
    #%
    #%   // big loop thru knot vector
    while b<m:
        #%   while (b < m)  {
        
    #%   }
    #% End big while loop                                      %   // end while loop
    #%
    #%   *nh = mh - ph - 1;
    #%
    #%   freevec2mat(ctrl);
    #%   freevec2mat(ictrl);
    #%   freematrix(bezalfs);
    #%   freematrix(bpts);
    #%   freematrix(ebpts);
    #%   freematrix(Nextbpts);
    #%   mxFree(alfs);
    #%
    #%   return(ierr);
    #% }
    return [ic, ik]
def bincoeff(n, k):

    # Local Variables: k, b, n
    # Function calls: bincoeff, factln, exp, floor
    #%  Computes the binomial coefficient.
    #%
    #%      ( n )      n!
    #%      (   ) = --------
    #%      ( k )   k!(n-k)!
    #%
    #%  b = bincoeff(n,k)
    #%
    #%  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
    #% double bincoeff(int n, int k)
    #% {
    b = np.floor((0.5+np.exp((factln(n)-factln(k)-factln((n-k))))))
    #%   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
    #% }
    return [b]
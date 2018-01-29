
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def basiskntins(deg, kc, kf):

    # Local Variables: kc, a, C, nk, alfa, kf, ss, nc, m, l, nf, n, ii, ik, r, u, jj, ind, b, nu, deg
    # Function calls: numel, speye, abs, new_knots, zeros, sparse, findspan, basiskntins
    #% Compute the coefficient matrix for non-uniform B-splines subdivision.
    #%
    #% This represents the B-spline basis given by a coarse knot vector
    #%  in terms of the B-spline basis of a finer knot vector.
    #%
    #% The function is implemented for the univariate case, based on
    #%  Algorithm A5.4 from 'The NURBS BOOK' pg164.
    #%
    #%
    #% Calling Sequence:
    #% 
    #%    S = basiskntins (deg, kc, kf);
    #%
    #%    INPUT:
    #%   
    #%      deg - degree of the first knot vector
    #%      kc  - coarse knot vector
    #%      kf  - fine knot vector
    #%   
    #%    OUTPUT:
    #%   
    #%      S - The matrix relating the two spaces, of size (deg-nu, deg-nt) 
    #%           with nu = numel(u)-deg-1, nt = numel(t)-deg-1
    #%   
    #%    Copyright (C) 2015, 2016 Rafael Vazquez
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
    nk = numel(kc)
    nc = nk-deg+1.
    u = new_knots(kc, kf)
    nu = numel(u)
    nf = nc+nu
    C = sparse(nf, nc)
    ik = np.zeros(1., (nk+nu))
    n = nc-1.
    r = nu-1.
    m = nc+deg
    a = findspan(n, deg, u[0], kc)
    b = findspan(n, deg, u[int(0)-1], kc)
    b = b+1.
    C[0:a-deg+1.,0:a-deg+1.] = speye((a-deg+1.))
    C[int(b+nu)-1:nc+nu,int(b)-1:nc] = speye((nc-b+1.))
    ik[0:a+1.] = kc[0:a+1.]
    ik[int(b+deg+nu+1.)-1:m+nu+1.] = kc[int(b+deg+1.)-1:m+1.]
    ii = b+deg-1.
    ss = ii+nu
    for jj in np.arange(r, (0.)+(-1.), -1.):
        ind = np.arange(a+1., (ii)+1)
        ind = ind[int((u[int((jj+1.))-1]<=kc[int((ind+1.))-1]))-1]
        C[int((ind+ss-ii-deg))-1,:] = 0.
        C[int((ind+ss-ii-deg))-1,int((ind-deg))-1] = speye(numel(ind))
        ik[int((ind+ss-ii+1.))-1] = kc[int((ind+1.))-1]
        ii = ii-numel(ind)
        ss = ss-numel(ind)
        C[int((ss-deg))-1,:] = C[int((ss-deg+1.))-1,:]
        for l in np.arange(1., (deg)+1):
            ind = ss-deg+l
            alfa = ik[int((ss+l+1.))-1]-u[int((jj+1.))-1]
            if np.abs(alfa) == 0.:
                C[int(ind)-1,:] = C[int((ind+1.))-1,:]
            else:
                alfa = matdiv(alfa, ik[int((ss+l+1.))-1]-kc[int((ii-deg+l+1.))-1])
                C[int(ind)-1,:] = np.dot(C[int(ind)-1,:], alfa)+np.dot(C[int((ind+1.))-1,:], 1.-alfa)
                
            
            
        ik[int((ss+1.))-1] = u[int((jj+1.))-1]
        ss = ss-1.
        
    return [C]
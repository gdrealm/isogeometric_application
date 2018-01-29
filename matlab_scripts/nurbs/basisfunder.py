
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def basisfunder(ii, pl, uu, u_knotl, nders):

    # Local Variables: j, ders, j1, j2, ii, right, saved, s2, s1, pk, pl, rk, nders, jj, u_knotl, a, d, temp, uu, k, dersv, i, r, u, ndu, left
    # Function calls: basisfunder, zeros, numel
    #% BASISFUNDER:  B-Spline Basis function derivatives.
    #%
    #% Calling Sequence:
    #% 
    #%   ders = basisfunder (ii, pl, uu, k, nd)
    #%
    #%    INPUT:
    #%   
    #%      ii  - knot span index (see findspan)
    #%      pl  - degree of curve
    #%      uu  - parametric points
    #%      k   - knot vector
    #%      nd  - number of derivatives to compute
    #%
    #%    OUTPUT:
    #%   
    #%      ders - ders(n, i, :) (i-1)-th derivative at n-th point
    #%   
    #%    Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72.
    #%
    #% See also: 
    #%
    #%    numbasisfun, basisfun, findspan
    #%
    #%    Copyright (C) 2009,2011 Rafael Vazquez
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
    dersv = np.zeros(numel(uu), (nders+1.), (pl+1.))
    for jj in np.arange(1., (numel(uu))+1):
        i = ii[int(jj)-1]+1.
        #%% convert to base-1 numbering of knot spans
        u = uu[int(jj)-1]
        ders = np.zeros((nders+1.), (pl+1.))
        ndu = np.zeros((pl+1.), (pl+1.))
        left = np.zeros((pl+1.))
        right = np.zeros((pl+1.))
        a = np.zeros(2., (pl+1.))
        ndu[0,0] = 1.
        for j in np.arange(1., (pl)+1):
            left[int((j+1.))-1] = u-u_knotl[int((i+1.-j))-1]
            right[int((j+1.))-1] = u_knotl[int((i+j))-1]-u
            saved = 0.
            for r in np.arange(0., (j-1.)+1):
                ndu[int((j+1.))-1,int((r+1.))-1] = right[int((r+2.))-1]+left[int((j-r+1.))-1]
                temp = matdiv(ndu[int((r+1.))-1,int(j)-1], ndu[int((j+1.))-1,int((r+1.))-1])
                ndu[int((r+1.))-1,int((j+1.))-1] = saved+np.dot(right[int((r+2.))-1], temp)
                saved = np.dot(left[int((j-r+1.))-1], temp)
                
            ndu[int((j+1.))-1,int((j+1.))-1] = saved
            
        for j in np.arange(0., (pl)+1):
            ders[0,int((j+1.))-1] = ndu[int((j+1.))-1,int((pl+1.))-1]
            
        for r in np.arange(0., (pl)+1):
            s1 = 0.
            s2 = 1.
            a[0,0] = 1.
            for k in np.arange(1., (nders)+1):
                #%compute kth derivative
                
            
        r = pl
        for k in np.arange(1., (nders)+1):
            for j in np.arange(0., (pl)+1):
                ders[int((k+1.))-1,int((j+1.))-1] = np.dot(ders[int((k+1.))-1,int((j+1.))-1], r)
                
            r = np.dot(r, pl-k)
            
        dersv[int(jj)-1,:,:] = ders
        
    #%!test
    #%! k    = [0 0 0 0 1 1 1 1];
    #%! p    = 3;
    #%! u    = rand (1);
    #%! i    = findspan (numel(k)-p-2, p, u, k);
    #%! ders = basisfunder (i, p, u, k, 1);
    #%! sumders = sum (squeeze(ders), 2);
    #%! assert (sumders(1), 1, 1e-15);
    #%! assert (sumders(2:end), 0, 1e-15);
    #%!test
    #%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
    #%! p    = 3;
    #%! u    = rand (1);
    #%! i    = findspan (numel(k)-p-2, p, u, k);
    #%! ders = basisfunder (i, p, u, k, 7); 
    #%! sumders = sum (squeeze(ders), 2);
    #%! assert (sumders(1), 1, 1e-15);
    #%! assert (sumders(2:end), zeros(rows(squeeze(ders))-1, 1), 1e-13);
    #%!test
    #%! k    = [0 0 0 0 1/3 2/3 1 1 1 1];
    #%! p    = 3;
    #%! u    = rand (100, 1);
    #%! i    = findspan (numel(k)-p-2, p, u, k);
    #%! ders = basisfunder (i, p, u, k, 7);
    #%! for ii=1:10
    #%!   sumders = sum (squeeze(ders(ii,:,:)), 2);
    #%!   assert (sumders(1), 1, 1e-15);
    #%!   assert (sumders(2:end), zeros(rows(squeeze(ders(ii,:,:)))-1, 1), 1e-13);
    #%! end
    #%! assert (ders(:, (p+2):end, :), zeros(numel(u), 8-p-1, p+1), 1e-13)
    #%! assert (all(all(ders(:, 1, :) <= 1)), true)
    return [dersv]
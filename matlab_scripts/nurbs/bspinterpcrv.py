
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspinterpcrv(Q, p, method):

    # Local Variables: A, ii, crv, d, n, Q, p, u, pnts, jj, x, y, span, z, knts, method
    # Function calls: basisfun, nrbmak, cumsum, sum, error, sqrt, nargin, ones, isempty, linspace, zeros, diff, size, bspinterpcrv, findspan, strcmpi
    #%
    #% BSPINTERPCRV: B-Spline interpolation of a 3d curve.
    #%
    #% Calling Sequence:
    #%
    #%   crv = bspinterpcrv (Q, p);
    #%   crv = bspinterpcrv (Q, p, method);
    #%   [crv, u] = bspinterpcrv (Q, p);
    #%   [crv, u] = bspinterpcrv (Q, p, method);
    #%
    #%    INPUT:
    #%
    #%      Q      - points to be interpolated in the form [x_coord; y_coord; z_coord].
    #%      p      - degree of the interpolating curve.
    #%      method - parametrization method. The available choices are:
    #%               'equally_spaced'
    #%               'chord_length'
    #%               'centripetal' (Default)
    #%
    #%    OUTPUT:
    #%
    #%      crv - the B-Spline curve.
    #%      u   - the parametric points corresponding to the interpolation ones.
    #%
    #%    See The NURBS book pag. 364 for more information.
    #%
    #%
    #% Copyright (C) 2015 Jacopo Corno
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
    #%
    if nargin<3. or isempty(method):
        method = 'centripetal'
    
    
    n = matcompat.size(Q, 2.)
    if strcmpi(method, 'equally_spaced'):
        u = np.linspace(0., 1., n)
    elif strcmpi(method, 'chord_length'):
        d = np.sum(np.sqrt(np.sum((np.diff(Q.conj().T).conj().T**2.), 1.)))
        u = np.zeros(1., n)
        u[1:n] = matdiv(np.cumsum(np.sqrt(np.sum((np.diff(Q, np.array([]), 2.)**2.), 1.))), d)
        #%    for ii = 2:n-1
        #%      u(ii) = u(ii-1) + norm (Q(:,ii) - Q(:,ii-1)) / d;
        #%    end
        u[int(0)-1] = 1.
        
    elif strcmpi(method, 'centripetal'):
        d = np.sum(np.sqrt(np.sqrt(np.sum((np.diff(Q.conj().T).conj().T**2.), 1.))))
        u = np.zeros(1., n)
        u[1:n] = matdiv(np.cumsum(np.sqrt(np.sqrt(np.sum((np.diff(Q, np.array([]), 2.)**2.), 1.)))), d)
        #%    for ii = 2:n-1
        #%      u(ii) = u(ii-1) + sqrt (norm (Q(:,ii) - Q(:,ii-1))) / d;
        #%    end
        u[int(0)-1] = 1.
        
    else:
        matcompat.error('BSPINTERPCRV: unrecognized parametrization method.')
        
    
    knts = np.zeros(1., (n+p+1.))
    for jj in np.arange(2., (n-p)+1):
        knts[int((jj+p))-1] = np.dot(1./p, np.sum(u[int(jj)-1:jj+p-1.]))
        
    knts[int(0-p)-1:] = np.ones(1., (p+1.))
    A = np.zeros(n, n)
    A[0,0] = 1.
    A[int(n)-1,int(n)-1] = 1.
    for ii in np.arange(2., (n-1.)+1):
        span = findspan(n, p, u[int(ii)-1], knts)
        A[int(ii)-1,int(span-p+1.)-1:span+1.] = basisfun(span, u[int(ii)-1], p, knts)
        
    x = linalg.solve(A, Q[0,:].conj().T)
    y = linalg.solve(A, Q[1,:].conj().T)
    z = linalg.solve(A, Q[2,:].conj().T)
    pnts = np.array(np.vstack((np.hstack((x.conj().T)), np.hstack((y.conj().T)), np.hstack((z.conj().T)), np.hstack((np.ones(matcompat.size(x.conj().T)))))))
    crv = nrbmak(pnts, knts)
    #%!demo
    #%! Q = [1 0 -1 -1 -2 -3;
    #%!      0 1  0 -1 -1 0;
    #%!      0 0  0  0  0 0];
    #%! p = 2;
    #%! crv = bspinterpcrv (Q, p);
    #%! 
    #%! plot (Q(1,:), Q(2,:), 'xk');
    #%! hold on; grid on;
    #%! nrbkntplot (crv);
    return [crv, u]
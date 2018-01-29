
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspinterpsurf(X, Y, Z, p, method):

    # Local Variables: A, ii, srf, d, P, span, m, n, Q, p, R, u, jj, v, Y, X, Z, knts, method
    # Function calls: squeeze, basisfun, nrbmak, sum, bspinterpsurf, sqrt, findspan, nargin, ones, zeros, linspace, diff, mean, isempty, size, cumsum, strcmp, strcmpi
    #%
    #% BSPINTERPSURF: B-Spline surface interpolation.
    #%
    #% Calling Sequence:
    #% 
    #%   srf = bspinterpsurf (Q, p, method);
    #%   
    #%    INPUT:
    #%   
    #%      X, Y, Z - grid of points to be interpolated. (See ndgrid)
    #%      p       - degree of the interpolating curve ([degree_x, degree_y]).
    #%      method  - parametrization method. The available choices are:
    #%                'equally_spaced'
    #%                'chord_length' (default)
    #%   
    #%    OUTPUT:
    #%   
    #%      srf - the B-Spline surface.
    #%   
    #%    See The NURBS book pag. 376 for more information. As of now only the
    #%    chord length method is implemented.
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
    if nargin<5. or isempty(method):
        method = 'chord_length'
    
    
    [n, m] = matcompat.size(X)
    Q = np.zeros(3., n, m)
    Q[0,:,:] = X
    Q[1,:,:] = Y
    Q[2,:,:] = Z
    if strcmpi(method, 'equally_spaced'):
        u = np.linspace(0., 1., n)
        v = np.linspace(0., 1., m)
    elif strcmp(method, 'chord_length'):
        u = np.zeros(m, n)
        for ii in np.arange(1., (m)+1):
            d = np.sum(np.sqrt(np.sum((np.diff(np.squeeze(Q[:,:,int(ii)-1]).conj().T).conj().T**2.), 1.)))
            u[int(ii)-1,1:n] = matdiv(np.cumsum(np.sqrt(np.sum((np.diff(Q[:,:,int(ii)-1], np.array([]), 2.)**2.), 1.))), d)
            #%      for jj = 2:n-1 
            #%        u(ii,jj) = u(ii,jj-1) + norm (Q(:,jj,ii) - Q(:,jj-1,ii)) / d;
            #%      end
            u[int(ii)-1,int(0)-1] = 1.
            
        u = np.mean(u)
        v = np.zeros(n, m)
        for ii in np.arange(1., (n)+1):
            d = np.sum(np.sqrt(np.sum((np.diff(np.squeeze(Q[:,int(ii)-1,:]).conj().T).conj().T**2.), 1.)))
            v[int(ii)-1,1:m] = matdiv(np.cumsum(np.sqrt(np.sum((np.diff(Q[:,int(ii)-1,:], np.array([]), 3.)**2.), 1.))), d)
            #%      for jj = 2:m-1 
            #%        v(ii,jj) = v(ii,jj-1) + norm (Q(:,ii,jj) - Q(:,ii,jj-1)) / d;
            #%      end
            v[int(ii)-1,int(0)-1] = 1.
            
        v = np.mean(v)
        
    
    #% TODO: implement centripetal method
    #% Compute knot vectors
    knts.cell[0] = np.zeros(1., (n+p[0]+1.))
    for jj in np.arange(2., (n-p[0])+1):
        knts.cell[0,int((jj+p[0]))-1][] = np.dot(1./p[0], np.sum(u[int(jj)-1:jj+p[0]-1.]))
        
    knts.cell[0,int(0-p[0])-1:][] = np.ones(1., (p[0]+1.))
    knts.cell[1] = np.zeros(1., (m+p[1]+1.))
    for jj in np.arange(2., (m-p[1])+1):
        knts.cell[1,int((jj+p[1]))-1][] = np.dot(1./p[1], np.sum(v[int(jj)-1:jj+p[1]-1.]))
        
    knts.cell[1,int(0-p[1])-1:][] = np.ones(1., (p[1]+1.))
    #% Interpolation
    R = np.zeros(matcompat.size(Q))
    P = np.zeros(4., n, m)
    for ii in np.arange(1., (m)+1):
        A = np.zeros(n, n)
        A[0,0] = 1.
        A[int(n)-1,int(n)-1] = 1.
        for jj in np.arange(2., (n-1.)+1):
            span = findspan(n, p[0], u[int(jj)-1], knts.cell[0])
            A[int(jj)-1,int(span-p[0]+1.)-1:span+1.] = basisfun(span, u[int(jj)-1], p[0], knts.cell[0])
            
        R[0,:,int(ii)-1] = linalg.solve(A, np.squeeze(Q[0,:,int(ii)-1]).conj().T)
        R[1,:,int(ii)-1] = linalg.solve(A, np.squeeze(Q[1,:,int(ii)-1]).conj().T)
        R[2,:,int(ii)-1] = linalg.solve(A, np.squeeze(Q[2,:,int(ii)-1]).conj().T)
        
    for ii in np.arange(1., (n)+1):
        A = np.zeros(m, m)
        A[0,0] = 1.
        A[int(m)-1,int(m)-1] = 1.
        for jj in np.arange(2., (m-1.)+1):
            span = findspan(m, p[1], v[int(jj)-1], knts.cell[1])
            A[int(jj)-1,int(span-p[1]+1.)-1:span+1.] = basisfun(span, v[int(jj)-1], p[1], knts.cell[1])
            
        P[0,int(ii)-1,:] = linalg.solve(A, np.squeeze(R[0,int(ii)-1,:]))
        P[1,int(ii)-1,:] = linalg.solve(A, np.squeeze(R[1,int(ii)-1,:]))
        P[2,int(ii)-1,:] = linalg.solve(A, np.squeeze(R[2,int(ii)-1,:]))
        
    P[3,:,:] = np.ones(n, m)
    #% Create B-Spline interpolant
    srf = nrbmak(P, knts)
    #%!demo
    #%! x = linspace (-3, 3, 40);
    #%! y = linspace (-3, 3, 40);
    #%! [X, Y] = meshgrid (x, y);
    #%! Z = peaks (X, Y);
    #%! 
    #%! srf1 = bspinterpsurf (X, Y, Z, [2 2], 'equally_spaced');
    #%! srf2 = bspinterpsurf (X, Y, Z, [2 2], 'chord_length');
    #%! figure
    #%! nrbkntplot(srf1)
    #%! title ('Approximation of the peaks functions, with the equally spaced method')
    #%! figure
    #%! nrbkntplot(srf2)
    #%! title ('Approximation of the peaks functions, with the chord length method')
    return [srf]
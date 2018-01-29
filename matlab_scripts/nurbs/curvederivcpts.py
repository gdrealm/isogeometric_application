
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def curvederivcpts(n, p, U, P, d, r1, r2):

    # Local Variables: tmp, d, r2, p, i, k, n, P, r, U, pk, r1
    # Function calls: nargin, curvederivcpts
    #% Compute control points of n-th derivatives of a B-spline curve.
    #% 
    #% usage: pk = curvederivcpts (n, p, U, P, d) 
    #%        pk = curvederivcpts (n, p, U, P, d, r1, r2) 
    #%
    #% If r1, r2 are not given, all the control points are computed.
    #%
    #%  INPUT:
    #%         n+1 = number of control points
    #%         p   = degree of the spline
    #%         d   = maximum derivative order (d<=p)
    #%         U   = knots
    #%         P   = control points
    #%         r1  = first control point to compute
    #%         r2  = auxiliary index for the last control point to compute
    #%  OUTPUT:
    #%         pk(k,i) = i-th control point of (k-1)-th derivative, r1 <= i <= r2-k
    #%
    #% Adaptation of algorithm A3.3 from the NURBS book, pg98.
    #%
    #%    Copyright (C) 2009 Carlo de Falco
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
    if nargin<=5.:
        r1 = 0.
        r2 = n
    
    
    r = r2-r1
    for i in np.arange(0., (r)+1):
        pk[0,int((i+1.))-1] = P[int((r1+i+1.))-1]
        
    for k in np.arange(1., (d)+1):
        tmp = p-k+1.
        for i in np.arange(0., (r-k)+1):
            pk[int((k+1.))-1,int((i+1.))-1] = matdiv(np.dot(tmp, pk[int(k)-1,int((i+2.))-1]-pk[int(k)-1,int((i+1.))-1]), U[int((r1+i+p+2.))-1]-U[int((r1+i+k+1.))-1])
            
        
    #%!test
    #%! line = nrbmak([0.0 1.5; 0.0 3.0],[0.0 0.0 1.0 1.0]);
    #%! pk   = curvederivcpts (line.number-1, line.order-1, line.knots,...
    #%!                        line.coefs(1,:), 2);
    #%! assert (pk, [0 3/2; 3/2 0], 100*eps);
    return [pk]
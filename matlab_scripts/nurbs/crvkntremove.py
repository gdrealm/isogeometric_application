
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def crvkntremove(crv, u, r, s, num, d):

    # Local Variables: Pw, crv, d, num, s, r, u, rcrv, U, t
    # Function calls: RemoveCurveKnot, crvkntremove, nrbmak
    #% 
    #% CRVKNTREMOVE: Remove one knot from the knot-vector of a NURBS curve.
    #% 
    #% Calling Sequence:
    #% 
    #%   [rcrv, remflag] = crvkntremove (crv, u, r, s, num, d);
    #% 
    #% INPUT:
    #% 
    #%   crv		: NURBS curve, see nrbmak.
    #% 
    #%   u           : knot to be removed.
    #% 
    #%   r           : index of the knot to be removed.
    #% 
    #%   s           : multiplicity of the knot to be removed.
    #% 
    #%   num         : number of knot removals requested.
    #%
    #%   d           : curve deviation tolerance.
    #%
    #% OUTPUT:
    #%
    #%   rcrv        : new NURBS structure for the curve with knot u remuved.
    #% 
    #%   t           : actual number of knot removals performed.
    #% 
    #% 
    #% 
    #% DESCRIPTION:
    #% 
    #%   Remove knot u from the NURBS curve crv at most num times. 
    #%   Check that the maximum deviation of the curve be less than d.
    #%   Based on algorithm A5.8 NURBS Book (pag183)
    #%  
    #% SEE ALSO:
    #% 
    #%   nrbkntins
    #%
    #%    Copyright (C) 2013 Jacopo Corno
    #%    Copyright (C) 2013 Carlo de Falco
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
    [U, Pw, t] = RemoveCurveKnot((crv.number), (crv.order-1.), (crv.knots), (crv.coefs), u, r, s, num, d)
    rcrv = nrbmak(Pw, U)
    #%!test 
    #%! crv  = nrbdegelev (nrbline (), 3);
    #%! acrv = nrbkntins (crv, [.11 .11 .11]);
    #%! [rcrv, t] = crvkntremove (acrv, .11, 8, 3, 3, 1e-10);
    #%! assert (crv.knots, rcrv.knots, 1e-10);
    #%! assert (t, 3);
    #%!test 
    #%! crv  = nrbcirc ();
    #%! acrv = nrbkntins (crv, [.3 .3]);
    #%! [rcrv, t] = crvkntremove (acrv, .3, 7, 2, 2, 1e-10);
    #%! assert (crv.knots, rcrv.knots, 1e-10);
    #%! assert (t, 2);
    return [rcrv, t]
def RemoveCurveKnot(n, p, U, Pw, u, r, s, num, d):

    # Local Variables: Pw, alfi, alfj, ii, num, remflag, fout, Pmax, U, jj, off, ord, d, last, temp, TOL, i, k, j, m, n, p, s, r, u, t, w, first
    # Function calls: min, max, sum, floor, sqrt, zeros, RemoveCurveKnot, norm, mod
    #% see algorithm A5.8 NURBS Book (pag183)
    w = matcompat.max(Pw[3,:])
    Pmax = matcompat.max(np.sqrt(np.sum((Pw**2.), 1.)))
    TOL = matdiv(np.dot(d, w), 1.+Pmax)
    m = n+p+1.
    ord = p+1.
    fout = (2.*r-s-p)/2.
    #% first control point out
    last = r-s
    first = r-p
    temp = np.zeros(4., (2.*p+1.))
    for t in np.arange(0., (num-1.)+1):
        off = first-1.
        #% diff in index between temp and P
        temp[:,0] = Pw[:,int(off)-1]
        temp[:,int((last+1.-off+1.))-1] = Pw[:,int((last+1.))-1]
        i = first
        j = last
        ii = 1.
        jj = last-off
        remflag = 0.
        while j-i > t:
            #% compute new control points for one removal step
            
        if j-i<=t:
            #% check if knot removable
        if linalg.norm((temp[:,int((ii-1.+1.))-1]-temp[:,int((jj+1.+1.))-1]))<=TOL:
            remflag = 1.
        else:
            alfi = matdiv(u-U[int(i)-1], U[int((i+ord+t))-1]-U[int(i)-1])
            if linalg.norm((Pw[:,int(i)-1]-alfi*temp[:,int((ii+t+1.+1.))-1]+(1.-alfi)*temp[:,int((ii-1.+1.))-1]))<=TOL:
                remflag = 1.
            
            
            #%if
            
        
        #%if  
        
        #%if
        if remflag == 0.:
            break
            #% cannot remove any more knots -> get out of for loop
        else:
            #% successful removal -> save new control points
            i = first
            j = last
            while j-i > t:
                Pw[:,int(i)-1] = temp[:,int((i-off+1.))-1]
                Pw[:,int(j)-1] = temp[:,int((j-off+1.))-1]
                i = i+1.
                j = j-1.
                
            
        
        #%if
        first = first-1.
        last = last+1.
        t = t+1.
        
    #% end of for loop
    if t == 0.:
        return [[rcrv, t]]
    
    
    #%if
    #% shift knots
    for k in np.arange(r+1., (m)+1):
        U[int((k-t))-1] = U[int(k)-1]
        
    U = U[0:0-t]
    j = np.floor(fout)
    i = j
    for k in np.arange(1., (t-1.)+1):
        if np.mod(k, 2.) == 1.:
            i = i+1.
        else:
            j = j-1.
            
        
        #%if
        
    #% shift points
    for k in np.arange(i+1., (n)+1):
        Pw[:,int(j)-1] = Pw[:,int(k)-1]
        j = j+1.
        
    Pw = Pw[:,0:0-t]
    return [[rcrv, t]]
    return [U, Pw, t]
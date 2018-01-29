
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% KNTREFINE: Refine a given knot vector by dividing each interval uniformly,
#%             maintaining the continuity in previously existing knots.
#%
#%   [rknots]                  = kntrefine (knots, n_sub, degree, regularity)
#%   [rknots, zeta]            = kntrefine (knots, n_sub, degree, regularity)
#%   [rknots, zeta, new_knots] = kntrefine (knots, n_sub, degree, regularity)
#%
#% INPUT:
#%
#%     knots:      initial knot vector.
#%     n_sub:      number of new knots to be added in each interval.
#%     degree:     polynomial degree of the refined knot vector
#%     regularity: maximum global regularity 
#%
#% OUTPUT:
#%
#%     rknots:    refined knot vector
#%     zeta:      refined knot vector without repetitions
#%     new_knots: new knots, to apply the knot insertion
#%
#% The regularity at the new inserted knots is the one given by the user.
#% At previously existing knots, the regularity is the minimum
#%  between the previous regularity, and the one given by the user.
#%  This ensures optimal convergence rates in the context of IGA.
#%
#% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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
def kntrefine(knots, n_sub, degree, regularity):

    # Local Variables: knots, old_mult, degree, insk, nz, varargout, regularity, aux_knots, n_sub, idim, new_knots, ik, min_mult, zeta, z, rknots, mult, deg
    # Function calls: repmat, max, sum, kntrefine, nargout, ones, linspace, numel, error, iscell, unique, vec
    if iscell(knots):
        if numel(n_sub) != numel(degree) or numel(n_sub) != numel(regularity) or numel(n_sub) != numel(knots):
            matcompat.error('kntrefine: n_sub, degree and regularity must have the same length as the number of knot vectors')
        
        
        aux_knots = knots
    else:
        if numel(n_sub) != numel(degree) or numel(n_sub) != numel(regularity) or numel(n_sub) != 1.:
            matcompat.error('kntrefine: n_sub, degree and regularity must have the same length as the number of knot vectors')
        
        
        aux_knots = cellarray(np.hstack((knots)))
        
    
    if nargout == 3.:
        for idim in np.arange(1., (numel(n_sub))+1):
            if degree[int(idim)-1]+1. != np.sum((aux_knots.cell[int(idim)-1] == aux_knots.cell[int(idim)-1,0]())):
                matcompat.error('kntrefine: new_knots is only computed when the degree is maintained')
            
            
            
        for idim in np.arange(1., (numel(n_sub))+1):
            min_mult = degree[int(idim)-1]-regularity[int(idim)-1]
            z = np.unique(aux_knots.cell[int(idim)-1])
            nz = numel(z)
            deg = np.sum((aux_knots.cell[int(idim)-1] == z[0]))-1.
            rknots.cell[int(idim)-1] = z[int(np.ones(1., (degree[int(idim)-1]+1.)))-1]
            new_knots.cell[int(idim)-1] = np.array([])
            for ik in np.arange(2., (nz)+1):
                insk = np.linspace(z[int((ik-1.))-1], z[int(ik)-1], (n_sub[int(idim)-1]+2.))
                insk = vec(matcompat.repmat(insk[1:0-1.], min_mult, 1.)).conj().T
                old_mult = np.sum((aux_knots.cell[int(idim)-1] == z[int(ik)-1]))
                mult = matcompat.max(min_mult, (degree[int(idim)-1]-deg+old_mult))
                rknots.cell[int(idim)-1] = np.array(np.hstack((rknots.cell[int(idim)-1], insk, z[int(np.dot(ik, np.ones(1., mult)))-1])))
                new_knots.cell[int(idim)-1] = np.array(np.hstack((new_knots.cell[int(idim)-1], insk, z[int(np.dot(ik, np.ones(1., (mult-old_mult))))-1])))
                
            zeta.cell[int(idim)-1] = np.unique(rknots.cell[int(idim)-1])
            
        if not iscell(knots):
            rknots = rknots.cell[0]
            zeta = zeta.cell[0]
            new_knots = new_knots.cell[0]
        
        
        varargout.cell[0] = rknots
        varargout.cell[1] = zeta
        varargout.cell[2] = new_knots
    else:
        for idim in np.arange(1., (numel(n_sub))+1):
            min_mult = degree[int(idim)-1]-regularity[int(idim)-1]
            z = np.unique(aux_knots.cell[int(idim)-1])
            nz = numel(z)
            deg = np.sum((aux_knots.cell[int(idim)-1] == z[0]))-1.
            rknots.cell[int(idim)-1] = z[int(np.ones(1., (degree[int(idim)-1]+1.)))-1]
            for ik in np.arange(2., (nz)+1):
                insk = np.linspace(z[int((ik-1.))-1], z[int(ik)-1], (n_sub[int(idim)-1]+2.))
                insk = vec(matcompat.repmat(insk[1:0-1.], min_mult, 1.)).conj().T
                old_mult = np.sum((aux_knots.cell[int(idim)-1] == z[int(ik)-1]))
                mult = matcompat.max(min_mult, (degree[int(idim)-1]-deg+old_mult))
                rknots.cell[int(idim)-1] = np.array(np.hstack((rknots.cell[int(idim)-1], insk, z[int(np.dot(ik, np.ones(1., mult)))-1])))
                
            zeta.cell[int(idim)-1] = np.unique(rknots.cell[int(idim)-1])
            
        if not iscell(knots):
            rknots = rknots.cell[0]
            zeta = zeta.cell[0]
        
        
        varargout.cell[0] = rknots
        if nargout == 2.:
            varargout.cell[1] = zeta
        
        
        
    
    return [varargout]
def vec(in):

    # Local Variables: in, v
    # Function calls: vec
    v = in.flatten(1)
    #%!shared nrbs
    #%!test
    #%! knots = {[0 0 1 1] [0 0 0 1 1 1]};
    #%! coefs(1,:,:) = [1 sqrt(2)/2 0; 2 sqrt(2) 0];
    #%! coefs(2,:,:) = [0 sqrt(2)/2 1; 0 sqrt(2) 2];
    #%! coefs(4,:,:) = [1 sqrt(2)/2 1; 1 sqrt(2)/2 1];
    #%! nrbs = nrbmak (coefs, knots);
    #%! nrbs = nrbkntins (nrbs, {[] [0.5 0.6 0.6]});
    #%! nrbs = nrbdegelev (nrbs, [0 1]);
    #%! nrbs = nrbkntins (nrbs, {[] [0.4]});
    #%! rknots = kntrefine (nrbs.knots, [1 1], [1 1], [0 0]);
    #%! assert (rknots{1} == [0 0 0.5 1 1]);
    #%! assert (rknots{2} == [0 0 0.2 0.4 0.45 0.5 0.55 0.6 0.8 1 1]);
    #%!
    #%!test
    #%! rknots = kntrefine (nrbs.knots, [1 1], [3 3], [0 0]);
    #%! assert (rknots{1}, [0 0 0 0 0.5 0.5 0.5 1 1 1 1]);
    #%! assert (rknots{2}, [0 0 0 0 0.2 0.2 0.2 0.4 0.4 0.4 0.45 0.45 0.45 0.5 0.5 0.5 0.55 0.55 0.55 0.6 0.6 0.6 0.8 0.8 0.8 1 1 1 1]);
    #%!
    #%!test
    #%! rknots = kntrefine (nrbs.knots, [1 1], [3 3], [2 2]);
    #%! assert (rknots{1}, [0 0 0 0 0.5 1 1 1 1]);
    #%! assert (rknots{2}, [0 0 0 0 0.2 0.4 0.45 0.5 0.5 0.55 0.6 0.6 0.6 0.8 1 1 1 1]);
    #%!
    #%!test
    #%! rknots = kntrefine (nrbs.knots, [1 1], [4 4], [0 0]);
    #%! assert (rknots{1}, [0 0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1 1]);
    #%! assert (rknots{2}, [0 0 0 0 0 0.2 0.2 0.2 0.2 0.4 0.4 0.4 0.4 0.45 0.45 0.45 0.45 0.5 0.5 0.5 0.5 0.55 0.55 0.55 0.55 0.6 0.6 0.6 0.6 0.8 0.8 0.8 0.8 1 1 1 1 1]);
    #%!
    #%!test
    #%! rknots = kntrefine (nrbs.knots, [1 1], [4 4], [3 3]);
    #%! assert (rknots{1}, [0 0 0 0 0 0.5 1 1 1 1 1]);
    #%! assert (rknots{2}, [0 0 0 0 0 0.2 0.4 0.4 0.45 0.5 0.5 0.5 0.55 0.6 0.6 0.6 0.6 0.8 1 1 1 1 1]);
    #%!
    #%!test
    #%! knots = [0 0 0 0 0.4 0.5 0.5 0.6 0.6 0.6 1 1 1 1];
    #%! rknots = kntrefine (knots, 1, 4, 3);
    #%! assert (rknots, [0 0 0 0 0 0.2 0.4 0.4 0.45 0.5 0.5 0.5 0.55 0.6 0.6 0.6 0.6 0.8 1 1 1 1 1]);
    return [v]
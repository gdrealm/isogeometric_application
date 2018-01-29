
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def aveknt(varargin):

    # Local Variables: ndim, nrb, onedim, idim, varargin, knt_aux, knt, pts, order
    # Function calls: false, nargin, reshape, sum, isfield, cell, aveknt, zeros, numel, error, iscell, repmat, true
    #% AVEKNT:  compute the knot averages (Greville points) of a knot vector
    #%
    #% Calling Sequence:
    #% 
    #%   pts = aveknt (knt, p)
    #%   pts = aveknt (nrb)
    #%   
    #%    INPUT:
    #%   
    #%      knt - knot sequence
    #%      p   - spline order (degree + 1)
    #%      nrb - NURBS structure (see nrbmak)
    #%   
    #%    OUTPUT:
    #%   
    #%      pts - average knots. If the input is a NURBS, it gives a cell-array,
    #%        with the average knots in each direction
    #%   
    #% See also: 
    #%
    #% Copyright (C) 2016 Rafael Vazquez
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
    if nargin == 1.:
        if isfield(varargin.cell[0], 'form'):
            nrb = varargin.cell[0]
            knt = nrb.knots
            order = nrb.order
        else:
            matcompat.error('The input should be a NURBS structure, or a knot vector and the order. See the help for details')
            
        
    elif nargin == 2.:
        knt = varargin.cell[0]
        order = varargin.cell[1]
        
    else:
        matcompat.error('The input should be a NURBS structure, or a knot vector and the order. See the help for details')
        
    
    onedim = false
    if not iscell(knt):
        knt = cellarray(np.hstack((knt)))
        onedim = true
    
    
    ndim = numel(knt)
    pts = cell(ndim, 1.)
    for idim in np.arange(1., (ndim)+1):
        if numel(knt.cell[int(idim)-1])<order[int(idim)-1]+1.:
            matcompat.error('The knot vector must contain at least p+2 knots, with p the degree')
        
        
        knt_aux = matcompat.repmat(knt.cell[int(idim)-1,1:0-1.](), 1., (order[int(idim)-1]-1.))
        knt_aux = np.array(np.vstack((np.hstack((knt_aux.flatten(1))), np.hstack((np.zeros((order[int(idim)-1]-1.), 1.))))))
        knt_aux = np.reshape(knt_aux, np.array([]), (order[int(idim)-1]-1.))
        pts.cell[int(idim)-1] = matdiv(np.sum(knt_aux.T, 1.), order[int(idim)-1]-1.)
        pts.cell[int(idim)-1] = pts.cell[int(idim)-1,0:0-order[int(idim)-1]+1.]()
        
    if onedim:
        pts = pts.cell[0]
    
    
    #%!test
    #%! knt = [0 0 0 0.5 1 1 1];
    #%! pts = aveknt (knt, 3);
    #%! assert (pts - [0 1/4 3/4 1] < 1e-14)
    #%!
    #%!test
    #%! knt = {[0 0 0 0.5 1 1 1] [0 0 0 0 1/3 2/3 1 1 1 1]};
    #%! pts = aveknt (knt, [3 4]);
    #%! assert (pts{1} - [0 1/4 3/4 1] < 1e-14);
    #%! assert (pts{2} - [0 1/9 1/3 2/3 8/9 1] < 1e-14);
    #%!
    #%!test
    #%! nrb = nrb4surf([0 0], [1 0], [0 1], [1 1]);
    #%! nrb = nrbkntins (nrbdegelev (nrb, [1 2]), {[1/2] [1/3 2/3]});
    #%! pts = aveknt (nrb);
    #%! assert (pts{1} - [0 1/4 3/4 1] < 1e-14);
    #%! assert (pts{2} - [0 1/9 1/3 2/3 8/9 1] < 1e-14);
    return [pts]
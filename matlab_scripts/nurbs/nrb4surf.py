
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def nrb4surf(p11, p12, p21, p22):

    # Local Variables: knots, srf, p11, p12, p21, p22, coefs
    # Function calls: nrb4surf, cat, nargin, length, ones, zeros, error
    #% 
    #% NRB4SURF: Constructs a NURBS bilinear surface.
    #% 
    #% Calling Sequence:
    #% 
    #%   srf = nrb4surf(p11,p12,p21,p22)
    #% 
    #% INPUT:
    #% 
    #%   p11		: Cartesian coordinate of the lhs bottom corner point.
    #% 
    #%   p12		: Cartesian coordinate of the rhs bottom corner point.
    #% 
    #%   p21		: Cartesian coordinate of the lhs top corner point.
    #%  
    #%   p22		: Cartesian coordinate of the rhs top corner point.
    #%
    #% OUTPUT:
    #% 
    #%   srf		: NURBS bilinear surface, see nrbmak. 
    #% 
    #% Description:
    #% 
    #%   Constructs a bilinear surface defined by four coordinates.
    #% 
    #%   The position of the corner points
    #% 
    #%          ^ V direction
    #%          |
    #%          ----------------
    #%          |p21        p22|
    #%          |              |
    #%          |    SRF       |
    #%          |              |
    #%          |p11        p12|
    #%          -------------------> U direction
    #% 
    #%
    #%    Copyright (C) 2000 Mark Spink
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
    #%!demo
    #%! srf = nrb4surf([0.0 0.0 0.5],[1.0 0.0 -0.5],[0.0 1.0 -0.5],[1.0 1.0 0.5]);
    #%! nrbplot(srf,[10,10]);
    #%! title('Construction of a bilinear surface.');
    #%! hold off
    return [srf]
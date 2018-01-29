
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% KNTBRKDEGMULT: Construct an open knot vector by giving the sequence of
#%                knots, the degree and the multiplicity.
#%
#%   knots = kntbrkdegreg (breaks, degree)
#%   knots = kntbrkdegreg (breaks, degree, mult)
#%
#% INPUT:
#%
#%     breaks:  sequence of knots.
#%     degree:  polynomial degree of the splines associated to the knot vector.
#%     mult:    multiplicity of the knots.
#%
#% OUTPUT:
#%
#%     knots:  knot vector.
#%
#% If MULT has as many entries as BREAKS, or as the number of interior
#%   knots, a different multiplicity will be assigned to each knot. If
#%   MULT is not present, it will be taken equal to 1.
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
def kntbrkdegmult(breaks, degree, mult):

    # Local Variables: knots, breaks, mult, degree
    # Function calls: false, do_kntbrkdegmult, kntbrkdegmult, nargin, num2cell, numel, error, cellfun, iscell
    if iscell(breaks):
        if nargin == 2.:
            mult = 1.
        
        
        if numel(breaks) != numel(degree) or numel(breaks) != numel(mult):
            matcompat.error('kntbrkdegmult: degree and multiplicity must have the same length as the number of knot vectors')
        
        
        degree = num2cell(degree)
        if not iscell(mult):
            mult = num2cell(mult)
        
        
        knots = cellfun(do_kntbrkdegmult, breaks, degree, mult, 'uniformoutput', false)
    else:
        if nargin == 2.:
            mult = 1.
        
        
        knots = do_kntbrkdegmult(breaks, degree, mult)
        
    
    return [knots]
def do_kntbrkdegmult(breaks, degree, mult):

    # Local Variables: knots, degree, mm, breaks, lm, sm, mults, mult
    # Function calls: sort, reshape, do_kntbrkdegmult, cumsum, sum, ones, zeros, numel, error, warning, any
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! knots = kntbrkdegmult (breaks, degree);
    #%! assert (knots, [0 0 0 0 1 2 3 4 4 4 4])
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! mult   = 2;
    #%! knots = kntbrkdegmult (breaks, degree, mult);
    #%! assert (knots, [0 0 0 0 1 1 2 2 3 3 4 4 4 4])
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! mult   = [1 2 3];
    #%! knots = kntbrkdegmult (breaks, degree, mult);
    #%! assert (knots, [0 0 0 0 1 2 2 3 3 3 4 4 4 4])
    #%!test
    #%! breaks = {[0 1 2 3 4] [0 1 2 3]};
    #%! degree = [3 2];
    #%! mult   = {[1 2 3] 2};
    #%! knots = kntbrkdegmult (breaks, degree, mult);
    #%! assert (knots, {[0 0 0 0 1 2 2 3 3 3 4 4 4 4] [0 0 0 1 1 2 2 3 3 3]})
    return [knots]

import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% KNTBRKDEGREG: Construct an open knot vector by giving the sequence of
#%                knots, the degree and the regularity.
#%
#%   knots = kntbrkdegreg (breaks, degree)
#%   knots = kntbrkdegreg (breaks, degree, regularity)
#%
#% INPUT:
#%
#%     breaks:     sequence of knots.
#%     degree:     polynomial degree of the splines associated to the knot vector.
#%     regularity: splines regularity.
#%
#% OUTPUT:
#%
#%     knots:  knot vector.
#%
#% If REGULARITY has as many entries as BREAKS, or as the number of interior
#%   knots, a different regularity will be assigned to each knot. If
#%   REGULARITY is not present, it will be taken equal to DEGREE-1.
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
def kntbrkdegreg(breaks, degree, reg):

    # Local Variables: knots, breaks, reg, degree
    # Function calls: false, do_kntbrkdegreg, nargin, num2cell, numel, error, cellfun, iscell, kntbrkdegreg
    if iscell(breaks):
        if nargin == 2.:
            reg = degree-1.
        
        
        if numel(breaks) != numel(degree) or numel(breaks) != numel(reg):
            matcompat.error('kntbrkdegreg: degree and regularity must have the same length as the number of knot vectors')
        
        
        degree = num2cell(degree)
        if not iscell(reg):
            reg = num2cell(reg)
        
        
        knots = cellfun(do_kntbrkdegreg, breaks, degree, reg, 'uniformoutput', false)
    else:
        if nargin == 2.:
            reg = degree-1.
        
        
        knots = do_kntbrkdegreg(breaks, degree, reg)
        
    
    return [knots]
def do_kntbrkdegreg(breaks, degree, reg):

    # Local Variables: knots, mults, breaks, reg, degree
    # Function calls: do_kntbrkdegreg, kntbrkdegmult, warning, numel, error, any
    if numel(breaks)<2.:
        matcompat.error('kntbrkdegreg: the knots sequence should contain at least two points')
    
    
    if numel(reg) == 1.:
        elif numel(reg) == numel(breaks):
            mults = degree-reg
            
    elif numel(reg) == numel(breaks)-2.:
        mults = np.array(np.hstack((-1., degree-reg, -1.)))
        
    else:
        matcompat.error('kntbrkdegreg: the length of mult should be equal to one or the number of knots')
        
    
    if np.any((reg<-1.)):
        matcompat.warning('kntbrkdegreg: for some knots the regularity is lower than -1')
    elif np.any((reg > degree-1.)):
        matcompat.error('kntbrkdegreg: the regularity should be lower than the degree')
        
    
    knots = kntbrkdegmult(breaks, degree, mults)
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! knots = kntbrkdegreg (breaks, degree);
    #%! assert (knots, [0 0 0 0 1 2 3 4 4 4 4])
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! reg    = 1;
    #%! knots = kntbrkdegreg (breaks, degree, reg);
    #%! assert (knots, [0 0 0 0 1 1 2 2 3 3 4 4 4 4])
    #%!test
    #%! breaks = [0 1 2 3 4];
    #%! degree = 3;
    #%! reg    = [0 1 2];
    #%! knots = kntbrkdegreg (breaks, degree, reg);
    #%! assert (knots, [0 0 0 0 1 1 1 2 2 3 4 4 4 4])
    #%!test
    #%! breaks = {[0 1 2 3 4] [0 1 2 3]};
    #%! degree = [3 2];
    #%! reg    = {[0 1 2] 0};
    #%! knots = kntbrkdegreg (breaks, degree, reg);
    #%! assert (knots, {[0 0 0 0 1 1 1 2 2 3 4 4 4 4] [0 0 0 1 1 2 2 3 3 3]})
    return [knots]

import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% KNTUNIFORM: generate uniform open knot vectors in the reference domain.
#%
#%   [csi, zeta] = kntuniform (num, degree, regularity)
#%
#% INPUT:
#%     
#%     num:        number of breaks (in each direction)
#%     degree:     polynomial degree (in each direction)
#%     regularity: global regularity (in each direction)
#%
#% OUTPUT:
#%
#%     csi:  knots
#%     zeta: breaks = knots without repetitions
#% 
#% Copyright (C) 2009, 2010 Carlo de Falco
#% Copyright (C) 2011 Rafael Vazquez
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
def kntuniform(num, degree, regularity):

    # Local Variables: csi, degree, rep, regularity, num, idim, zeta
    # Function calls: linspace, zeros, kntuniform, numel, error
    if numel(num) != numel(degree) or numel(num) != numel(regularity):
        matcompat.error('kntuniform: num, degree and regularity must have the same length')
    else:
        for idim in np.arange(1., (numel(num))+1):
            zeta.cell[int(idim)-1] = np.linspace(0., 1., num[int(idim)-1])
            rep = degree[int(idim)-1]-regularity[int(idim)-1]
            if rep > 0.:
                else:
                    matcompat.error('kntuniform: regularity requested is too high')
                    
            
            
            
        if numel(num) == 1.:
            csi = csi.cell[0]
            zeta = zeta.cell[0]
        
        
        
    
    return [csi, zeta]
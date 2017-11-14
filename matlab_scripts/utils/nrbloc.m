function t = nrbloc(crv,d,x,params)
%
% NRBLOC:

% INPUT:
%
%   crv: NURBS structure of the curve
%   d: dimension to search for local point
%   x: global coordinate of point to search for local coordinate
%
% OUTPUT:
%   t: the local coordinate of x in the curve

% default parameters
if ~isfield(params,'max_it')
    params.max_it = 30;
end
if ~isfield(params,'tol')
    params.tol = 1.0e-10;
end
converged = 0;

% choose an initial value
t = 0.0;
dcrv = nrbderiv(crv);

% compute the initial residual
[p, dp] = nrbdeval(crv,dcrv,t);
res = x - p(d);

% iteration to find local point
it = 0;
while(converged==0 && it<params.max_it)
    if(res < params.tol)
        converged = 1;
        break;
    end

    t = t + res / dp(d);

    [p, dp] = nrbdeval(crv,dcrv,t);
    res = x - p(d);

    it = it + 1;
end

if(it==params.max_it && converged==0)
    error(['The iteration does not converge in' params.max_it ' steps'])
end

end


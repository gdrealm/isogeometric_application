function t = nrbloc1(knots,p,w,v,x,params)
%
% NRBLOC1: 1-dimension version of nrbloc

% INPUT:
%
%   knots: knot vector
%   p: precise order of the nurbs curve
%   w: vector of control point weight
%   v: vector of control point value (must be already weighted)
%   x: global point want to find local coordinate from
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
coefs = [v;w];
[dcoefs,dknots] = bspderiv(p,coefs,knots);

% compute the initial residual
cp = bspeval(p,coefs,knots,t);
pnt = cp(1)/cp(2);
res = x - pnt;

% iteration to find local point
it = 0;
while(converged==0 && it<params.max_it)
    if(res < params.tol)
        converged = 1;
        break;
    end

    dcp = bspeval(p-1,dcoefs,dknots,t);
    dpnt = (dcp(1)-dcp(2)*pnt)/cp(2);

    t = t + res / dpnt;

    cp = bspeval(p,coefs,knots,t);
    pnt = cp(1)/cp(2);
    res = x - pnt;

    it = it + 1;
end

if(it==params.max_it && converged==0)
    error('The iteration does not converge in 30 steps')
end

end


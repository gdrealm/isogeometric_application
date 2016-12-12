function [C,nt,Ubar] = bezierExtractionLocal1DCell(Xi,Uxi,p)
% Compute the Bezier extraction of T-splines basis function over a cell; this cell must lie in one knot span of local knot vector
% Input:
%   Xi      local knot vector
%   Uxi     cell (a,b)
%   p       curve degree
%
% Output:
%   C       matrix of size (p+1) contain the rows of bezier extraction of the T-splines basis function over the cell
%   nt      span of the cell in Ubar
%   Ubar    knot vector including the cell

% do some check
if length(Xi) ~= p+2
    error('local knot vector must be of length p + 2')
end

if length(Uxi) ~= 2
    error('1D cell must contain 2 end')
end

if Uxi(1) >= Uxi(2)
    error('cell must be nonzero positive length')
end

% locate the span of the cell
span = 1;
while Uxi(1) >= Xi(span)
    span = span + 1;
end
if Uxi(2) > Xi(span)
    error('cell must lie in one knot span')
end
span = span - 1;

% compute the inserted knots 
ins_knots = [];
if Uxi(1) ~= Xi(span)
    ins_knots = [ins_knots Uxi(1)];
end
if Uxi(2) ~= Xi(span+1)
    ins_knots = [ins_knots Uxi(2)];
end

% compute the Bezier extraction operator
[D,nb,Ubar] = bezierExtractionLocal1D(Xi,ins_knots,p);
D
Ubar = unique(Ubar);
for i = 1:length(Ubar)-1
    if Uxi(1)>=Ubar(i) && Uxi(2)<=Ubar(i+1)
        nt = i;
        break
    end
end
C = D(nt,:);


%% compute the full inserted knot vector for Bezier decomposition
% Input
%       U    inserted knots
%       p    degree of the basis function
% Output
%      Ub      full inserted knots
function Ub = compute_bezier_inserted_knots(U,p)
% create the vector of unique inserted knots
Un = unique(U);

% replicate each knots p times
Ub = zeros(1,p*length(Un));
for i = 1:length(Un)
    Ub((i-1)*p+1:i*p) = Un(i);
end


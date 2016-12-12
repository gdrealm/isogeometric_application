%% compute the extended knot vector given the local knot vector, assuming the local knot vector is in ascending order
% Input
%      Xi    local knot vector
%       p    degree of the basis function
% Output
%    Ubar    extended knot vector
%      nt    index of this basis function w.r.t the extended knot vector
function [Ubar,nt] = compute_extended_knot_vector(Xi,p)
% count the multiplicity of the first knot
n = length(Xi);
a = 0;
for i = 1:n
    if Xi(i) == Xi(1)
        a = a + 1;
    else
        break;
    end
end
nt = p - a + 1;

% count the multiplicity of the last knot
b = 0;
for i = n:-1:1
    if Xi(i) == Xi(n)
        b = b + 1;
    else
        break;
    end
end

% compute the extended knot vector
Ubar = [Xi(1)*ones(1,nt) Xi Xi(n)*ones(1,p-b+1)];

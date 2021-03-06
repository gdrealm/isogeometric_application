function [C,nb,Ubar] = bezierExtractionLocal1D(Xi,U,p)
% Input:
%   Xi      local knot vector
%   U       inserted knots
%   p       curve degree
%
% Output:
%   C       matrix of size nb x (p+1) contain the rows of bezier extraction of the T-splines basis function at each knot span of the filled extended knot vector
%   nb      number of spans of the filled extended knot vector
%   Ubar    filled extended knot vector

% do some check
if length(Xi) ~= p+2
    error('local knot vector must be of length p + 2')
end

% compute the extended knot vector
[Ubar,nt] = compute_extended_knot_vector(Xi,p);

% count the multiplicity of inner knots
i = p+2;
Ud = [];
Um = [];
num_inner_knots = 0;
while i <= length(Ubar)-p-1
    mult = 0;
    while Ubar(i+mult) == Ubar(i)
        mult = mult + 1;
    end
    num_inner_knots = num_inner_knots + 1;
    Ud(num_inner_knots) = Ubar(i);
    Um(num_inner_knots) = mult;
    i = i + mult;
end

% compute the inserted knots
ins_knots = [];
cnt = 1;
for i = 1:num_inner_knots
    if Um(i) < p
        ins_knots(cnt:cnt+(p-Um(i)-1)) = Ud(i)
        cnt = cnt + (p-Um(i));
    end
end

% append the external knots
for i = 1:length(U)
    found = 0;
    for j = 1:length(Ud)
        if Ud(j) == U(i)
            found = 1;
            break
        end
    end
    if found == 0
        num_inner_knots = num_inner_knots + 1;
        ins_knots(length(ins_knots)+1:length(ins_knots)+p) = U(i);
    end
end

% compute the full extraction operator
D = eye(length(Ubar)-p-1);
for i = 1:length(ins_knots)
    [Di,Ubar] = bspkntins_coeffs_1(p,Ubar,ins_knots(i));
    D = D*Di;
end

% extract the local extraction operator
Tmp = D(nt+1,:);
nb = num_inner_knots+1;
j = 1;
C = zeros(nb,p+1);
for i = 1:nb
    C(i,:) = Tmp(j:j+p);
    j = j + p;
end


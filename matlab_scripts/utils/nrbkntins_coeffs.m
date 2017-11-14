%% compute the matrix of linear combination coefficients of NURBS refinement
function [D,new_knot_vector,new_w_vector] = nrbkntins_coeffs(nurbs,iknots)
% firstly compute the refinement matrix for B-spline refinement
[D,new_knot_vector] = bspkntins_coeffs(nurbs,iknots);
% D

% secondly compute the new weight vector
w_vector = nurbs.coefs(4,:);
new_w_vector = D'*w_vector';

% thirdly compute the coefficient matrix
for i = 1:size(D,1)
    for j = 1:size(D,2)
        D(i,j) = D(i,j)*w_vector(i)/new_w_vector(j);
    end
end

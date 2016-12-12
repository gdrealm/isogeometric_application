%% compute the refinement matrix for local knot vector
function [D,new_knots] = bspkntins_coeffs_2(knots,iknots,p)
[full_knot_vector,nt] = compute_extended_knot_vector(knots,p);
% full_knot_vector
% nt
n = length(full_knot_vector)-p-1;
D = eye(n);
for i = 1:length(iknots)
    [Di,full_knot_vector] = bspkntins_coeffs_1(p,full_knot_vector,iknots(i));
    D = D*Di;
end
% full_knot_vector
% D(nt+1,:)
D = D(nt+1,nt+1:nt+1+length(iknots));
new_knots = full_knot_vector(nt+1:nt+1+length(iknots)+p+1);

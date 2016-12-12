function [C,new_knot_vector] = bspkntins_coeffs_1(p,knot_vector,k)
% compute number of basis function
n = length(knot_vector)-p-1;

% find the span of the knot
s = findspan(n-1,p,k,knot_vector)+1;

% form the new knot
new_knot_vector= [knot_vector(1:s) k knot_vector(s+1:end)];

% initialize and compute matrix C
C = zeros(n,n+1);
for i = 1:s-p-1
    C(i,i) = 1.0;
end
for i = s-p:s
    C(i,i) = (k-new_knot_vector(i))/(new_knot_vector(i+p+1)-new_knot_vector(i));
    C(i,i+1) = (new_knot_vector(i+p+2)-k)/(new_knot_vector(i+p+2)-new_knot_vector(i+1));
end
for i = s+1:n
    C(i,i+1) = 1.0;
end

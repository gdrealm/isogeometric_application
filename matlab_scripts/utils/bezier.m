function y = bezier(u,i,n)
% evaluate Bezier shape function
% Input
%   u    evaluation point
%   i    index of the shape function, 0 <= i <= n
%   n    degree of the shape function
y = nchoosek(n,i) * u^i * (1-u)^(n-i);

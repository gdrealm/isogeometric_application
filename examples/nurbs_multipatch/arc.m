knots = [ 0 0 0 1 1 1];
coefs = zeros(4,3);
coefs(:,1) = [ 1 0 0 1];
coefs(:,2) = [ 0.866025 3.06162e-17 0.5 0.866025];
coefs(:,3) = [ 0.5 5.30288e-17 0.866025 1];
nurbs = nrbmak(coefs,knots);


knots = {};
knots{1} = [ 0 0 0 1 1 1];
knots{2} = [ 0 0 0 0.5 1 1 1];
coefs = zeros(4,3,4);
coefs(:,1,1) = [ 0 0 0 1];
coefs(:,2,1) = [ 1.5 0 0 1];
coefs(:,3,1) = [ 3 0 0 1];
coefs(:,1,2) = [ 0 4 0 1];
coefs(:,2,2) = [ 1.5 0.75 0 1];
coefs(:,3,2) = [ 3 0.41421 0 1];
coefs(:,1,3) = [ 0 4 0 1];
coefs(:,2,3) = [ 3.25 2.5 0 1];
coefs(:,3,3) = [ 3.58578 1 0 1];
coefs(:,1,4) = [ 4 4 0 1];
coefs(:,2,4) = [ 4 2.5 0 1];
coefs(:,3,4) = [ 4 1 0 1];
nurbs = nrbmak(coefs,knots);


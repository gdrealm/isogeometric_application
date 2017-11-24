knots = {};
knots{1} = [ 0 0 0 0 1 1 1 1];
knots{2} = [ 0 0 1 1];
coefs = zeros(4,4,2);
coefs(:,1,1) = [ 0 0 0 1];
coefs(:,2,1) = [ 0.333333 0 0 1];
coefs(:,3,1) = [ 0.666667 0 0 1];
coefs(:,4,1) = [ 1 0 0 1];
coefs(:,1,2) = [ 0 1 0 1];
coefs(:,2,2) = [ 0.333333 1 0 1];
coefs(:,3,2) = [ 0.666667 1 0 1];
coefs(:,4,2) = [ 1 1 0 1];
nurbs = nrbmak(coefs,knots);


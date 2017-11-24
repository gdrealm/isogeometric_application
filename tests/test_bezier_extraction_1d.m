%% generate a curve and its bezier extraction operator
knots = [0.00  0.00  0.00  0.00  0.25  0.50  0.75  1.00  1.00  1.00  1.00];
coefs = [
 0.00  0.00  1.00  2.00  3.00  4.00  5.00
 0.00  1.00  1.00  0.00  2.00  2.25  1.00
 0.00  0.00  0.00  0.00  0.00  0.00  0.00
 1.00  1.00  1.00  1.00  1.00  1.00  1.00];
 %1.00  0.90  0.80  0.70  0.60  0.50  0.40
 
nurbs = nrbmak(coefs, knots);

[C nb] = bezierExtraction(nurbs.knots, nurbs.order - 1);
nb
C



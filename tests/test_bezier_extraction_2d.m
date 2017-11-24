geom = geo_load('geo_ring.txt');
nurbs = geom.nurbs;

nurbs = p_refine(nurbs, [2 2]);

%nurbs.coefs

nurbs = h_refine(nurbs, [2 2], [1 1]);

% nrbplot(geom.nurbs, [100 100]);

nurbs.knots{1}
nurbs.knots{2}

%% export to layer file
[C, Cet, Cxi] = bezierExtraction2D(nurbs.knots{1}, ...
    nurbs.knots{2}, nurbs.order(1) - 1, ...
    nurbs.order(2) - 1);
C
Cet
Cxi


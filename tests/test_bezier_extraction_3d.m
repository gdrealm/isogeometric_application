geom = geo_load('slope3d.txt');
nurbs = geom.nurbs;

nurbs.knots{1}
nurbs.knots{2}
nurbs.knots{3}

%% export to layer file
[C, Cet, Cxi, Cze] = bezierExtraction3D(nurbs.knots{1}, ...
    nurbs.knots{2}, nurbs.knots{3}, nurbs.order(1) - 1, ...
    nurbs.order(2) - 1, nurbs.order(3) - 1);
C

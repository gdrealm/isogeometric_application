geom = geo_load('infinite_plate.txt');
nurbs = geom.nurbs;

% geom2 = p_refine(geom1, [2 2]);
% 
% geom2.nurbs.coefs
% 
% geom = h_refine(geom2, [2 2], [1 1]);

% nrbplot(nurbs, [100 100]);
%params.axis = 'off';
%plot_ctrl_points_2d(nurbs,params);

nurbs.knots{1}
nurbs.knots{2}

[C, Cet, Cxi] = bezierExtraction2D(nurbs.knots{1}, ...
    nurbs.knots{2}, nurbs.order(1) - 1, ...
    nurbs.order(2) - 1);

C


addpath('/home/hbui/kratos_bcn2/applications/isogeometric_application/matlab_scripts/nurbs');
addpath('/home/hbui/kratos_bcn2/applications/isogeometric_application/matlab_scripts/geopde/geometry');
addpath('/home/hbui/kratos_bcn2/applications/isogeometric_application/matlab_scripts/utils');

geom1 = geo_load('patch1.txt');
nurbs1 = geom1.nurbs;
params.label = 'off';
%plot_ctrl_points_2d(nurbs1,params);
nurbs1_1 = nrbkntins(nurbs1, {[0.5] [0.5]});

geom2 = geo_load('patch2.txt');
nurbs2 = geom2.nurbs;
nurbs2_1 = nrbkntins(nurbs2, {[] [0.5]});



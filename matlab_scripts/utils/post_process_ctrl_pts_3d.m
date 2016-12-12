%% post process the control point results to the isogeometric mesh
% this only works with single patch
%
%
%
function post_process_ctrl_pts_3d(u,geom,method_data,vtk_pts,fn,dofname)

% Construct msh structure
rule     = msh_gauss_nodes(method_data.nquad);
[qn, qw] = msh_set_quad_nodes(geom.nurbs.knots,rule);
msh      = msh_3d(geom.nurbs.knots,qn,qw,geom);

% Construct space structure
sp_scalar = sp_nurbs_3d(geom.nurbs,msh);
sp = sp_vector_3d(sp_scalar,sp_scalar,sp_scalar,msh);
sp.ndof
%clear sp_scalar

%sp2 = sp_vector_2d(sp_scalar,sp_scalar,msh);
%sp2.ndof

sp_to_vtk(u,sp,geom,vtk_pts,sprintf('%s_%s',fn,dofname),dofname);

end


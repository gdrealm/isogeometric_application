clc
close all

%% add path
addpath('../../nurbs-util');
addpath('../../utils');
addpath('../../mesh_data_set_2--half_straight_tunnel');

geometry = geo_load('tunnel_virgin2.txt');

% Material properties
E  = 5.0e+7;
nu = 0.3; 
lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Gravity
rho = 2400.0;
g = -9.81;
fx = @(x, y, z) zeros(size(x));
fy = @(x, y, z) zeros(size(x));
fz = @(x, y, z) rho * g * ones(size(x));
f = @(x, y, z) cat(1, ...
                    reshape (fx (x,y,z), [1, size(x)]), ...
                    reshape (fy (x,y,z), [1, size(x)]), ...
                    reshape (fz (x,y,z), [1, size(x)]));
                
% Number of points for Gaussian quadrature rule
nquad = [3 3 3];

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_3d (geometry.nurbs.knots, qn, qw, geometry);

% Construct space structure
sp_scalar = sp_nurbs_3d (nurbs, msh);
sp = sp_vector_3d (sp_scalar, sp_scalar, sp_scalar, msh);
clear sp_scalar

% Assemble the matrices
mat    = op_su_ev_tp (sp, sp, msh, lambda_lame, mu_lame); 
rhs    = op_f_v_tp (sp, msh, f);

% Dirichlet boundary condition
drchlt_dofs = [1:12 13:24 25 26 31 32];

% Solve
int_dofs = setdiff (1:sp.ndof, [drchlt_dofs]);
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
mat(int_dofs, int_dofs)
rhs(int_dofs)
u(int_dofs)'

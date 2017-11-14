function write_layer_bezier_2d_boundary_marker(fid, name, nurbs)
%
% write the Bezier boundary
%
% Input:
%   fid     file handler
%   name    layer name
%   nurbs   nurbs geometry

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

%%
fprintf(fid, '\t\t# layer boundary marker\n');
fprintf(fid, '\t\tboundary_marker = {}\n');

% extract the boundary control point indices
boundary_uv = zeros(u_dim,v_dim);
boundary_u0 = zeros(1,v_dim);
boundary_u1 = zeros(1,v_dim);
boundary_v0 = zeros(1,u_dim);
boundary_v1 = zeros(1,u_dim);
cnt = 1;
for j = 1:v_dim
    for k = 1:u_dim
        boundary_uv(k,j) = cnt;

        if k==1
            boundary_u0(j) = cnt;
        elseif k==u_dim
            boundary_u1(j) = cnt;
        end

        if j==1
            boundary_v0(k) = cnt;
        elseif j==v_dim
            boundary_v1(k) = cnt;
        end

        cnt = cnt + 1;
    end
end

fprintf(fid, '\t\tboundary_marker[''uv''] = ');
print_integer_matrix(fid,boundary_uv);

fprintf(fid, '\t\tboundary_marker[''u0''] = ');
print_integer_matrix(fid,boundary_u0);
fprintf(fid, '\t\tboundary_marker[''u1''] = ');
print_integer_matrix(fid,boundary_u1);
fprintf(fid, '\t\tboundary_marker[''v0''] = ');
print_integer_matrix(fid,boundary_v0);
fprintf(fid, '\t\tboundary_marker[''v1''] = ');
print_integer_matrix(fid,boundary_v1);

fprintf(fid, '\t\tself.layer_boundary_marker[''%s''] = boundary_marker\n\n', name);


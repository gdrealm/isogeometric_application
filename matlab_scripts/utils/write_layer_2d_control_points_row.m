function write_layer_bezier_2d_control_points_row(fid, name, nurbs, row_name)
%
% write the control points row for NURBS in 2D
% each row will be written line by line to list
%
% Input:
%   fid     file handler
%   name    layer name
%   nurbs   nurbs geometry
%   row_name   'u' or 'v'

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

%%
fprintf(fid, '\t\t# control points row\n');
fprintf(fid, '\t\tcontrol_points_row_%s = []\n', row_name);

% extract the control points row
if row_name == 'u'
    for j = 1:v_dim
        fprintf(fid, '\t\trow = [');
        cnt = 1 + (j-1)*u_dim;
        for k = 1:u_dim
            fprintf(fid, '%d, ', cnt);
            cnt = cnt + 1;
        end
        fprintf(fid, ']\n');
        fprintf(fid, '\t\tcontrol_points_row_%s.append(row)\n', row_name);
    end
elseif row_name == 'v'
    for k = 1:u_dim
        fprintf(fid, '\t\trow = [');
        cnt = k;
        for j = 1:v_dim
            fprintf(fid, '%d, ', cnt);
            cnt = cnt + u_dim;
        end
        fprintf(fid, ']\n');
        fprintf(fid, '\t\tcontrol_points_row_%s.append(row)\n', row_name);
    end
end

fprintf(fid, '\t\tself.control_points_row_%s[''%s''] = control_points_row_%s\n\n', row_name, name, row_name);


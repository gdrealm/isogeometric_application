function plot_ctrl_points_3d(nurbs,params)
% Remarks: it works nicely in the script, but in the command line may exhibit weird behaviour
% n = size(x(:), 1);
% hold on
% for i = 1:n
%     text(x(i), y(i), z(i), num2str(i));
% end

axis equal;
hold on;

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);
w_dim = sizes(4);

u_color = 'black';
v_color = 'cyan';
w_color = 'red';

if ~isfield(params,'point_color')
    params.point_color = 'blue';
end

if ~isfield(params,'label')
    params.label = 'off';
end

if ~isfield(params,'axis')
    params.label = 'off';
end

if ~isfield(params,'text_dc')
    params.text_dc = 1.01;
end

%%
cnt = 1;
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            scatter3(point(1), point(2), point(3), params.point_color);
            if strcmp(params.label,'on')
                text(point(1)*params.text_dc, point(2)*params.text_dc, point(3)*params.text_dc, num2str(cnt));
            end
            cnt = cnt + 1;
            if k > 1
                if k == 2
                    L = quiver3(old_point_u(1),old_point_u(2),old_point_u(3),point(1)-old_point_u(1),point(2)-old_point_u(2),point(3)-old_point_u(3));
%                    set(L, 'markersize', 10.0);
                else
                    L = line([old_point_u(1) point(1)], [old_point_u(2) point(2)], [old_point_u(3) point(3)]);
                end
                set(L, 'color', u_color);
            end
            old_point_u = point;
        end
    end
end

%%
for i = 1:w_dim
    for k = 1:u_dim
        for j = 1:v_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            if j > 1
                if j == 2
                    L = quiver3(old_point_v(1),old_point_v(2),old_point_v(3),point(1)-old_point_v(1),point(2)-old_point_v(2),point(3)-old_point_v(3));
                else
                    L = line([old_point_v(1) point(1)], [old_point_v(2) point(2)], [old_point_v(3) point(3)]);
                end
                set(L, 'color', v_color);
            end
            old_point_v = point;
        end
    end
end

%%
for j = 1:v_dim
    for k = 1:u_dim
        for i = 1:w_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            if i > 1
                if i == 2
                    L = quiver3(old_point_w(1),old_point_w(2),old_point_w(3),point(1)-old_point_w(1),point(2)-old_point_w(2),point(3)-old_point_w(3));
                else
                    L = line([old_point_w(1) point(1)], [old_point_w(2) point(2)], [old_point_w(3) point(3)]);
                end
                set(L, 'color', w_color);
            end
            old_point_w = point;
        end
    end
end

u = plot(0,0,u_color);
v = plot(0,0,v_color);
w = plot(0,0,w_color);
legend([u,v,w], 'u-dim', 'v-dim', 'w-dim');

xlabel('x');
ylabel('y');
zlabel('z');

if strcmp(params.axis,'on')
    Lx = line([0 1], [0 0], [0 0]);
    text(1, 0, 0, 'X');
    set(Lx, 'color', 'magenta');

    Ly = line([0 0], [0 1], [0 0]);
    text(0, 1, 0, 'Y');
    set(Ly, 'color', 'magenta');

    Lz = line([0 0], [0 0], [0 1]);
    text(0, 0, 1, 'Z');
    set(Lz, 'color', 'magenta');
end


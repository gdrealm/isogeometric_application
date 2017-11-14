function plot_ctrl_points_2d(nurbs,params)

axis equal
hold on

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

u_color = 'black';
v_color = 'cyan';

if ~isfield(params,'point_color')
    params.point_color = 'blue';
end

if ~isfield(params,'label')
    params.label = 'off';
end

%%
cnt = 1;
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j);
            point(1:3) = point(1:3) / point(4);
            scatter3(point(1), point(2), point(3), params.point_color);
            if strcmp(params.label,'on')
                text(point(1), point(2), point(3), num2str(cnt), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'cap');
            end
            cnt = cnt + 1;
            if k > 1
                if k == 2
                    L = quiver3(old_point_u(1), old_point_u(2), old_point_u(3), point(1)-old_point_u(1), point(2)-old_point_u(2), point(3)-old_point_u(3));
                else
                    L = line([old_point_u(1) point(1)], [old_point_u(2) point(2)], [old_point_u(3) point(3)]);
                end
                set(L, 'color', u_color);
            end
            old_point_u = point;
        end
    end


%%
    for k = 1:u_dim
        for j = 1:v_dim
            point = nurbs.coefs(:, k, j);
            point(1:3) = point(1:3) / point(4);
            if j > 1
                if j == 2
                    L = quiver3(old_point_v(1), old_point_v(2), old_point_v(3), point(1)-old_point_v(1), point(2)-old_point_v(2), point(3)-old_point_v(3));
                else
                    L = line([old_point_v(1) point(1)], [old_point_v(2) point(2)], [old_point_v(3) point(3)]);
                end
                set(L, 'color', v_color);
            end
            old_point_v = point;
        end
    end

u = plot(0,0,u_color);
v = plot(0,0,v_color);
legend([u,v], 'u-dim', 'v-dim');

if ~isfield(params,'axis')
    params.axis = 'off';
end

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


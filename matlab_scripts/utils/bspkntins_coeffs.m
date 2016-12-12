%% compute the matrix of linear combination coefficients of b-splines refinement in 1D
% Ref(partially): http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS-knot-insert.html
function [D,knot_vector] = bspkntins_coeffs(nurbs,iknots)
knot_vector = nurbs.knots;
if ~iscell(knot_vector) % curve
    D = eye(nurbs.number);
    for i = 1:length(iknots)
        [Di,knot_vector] = bspkntins_coeffs_1(nurbs.order-1,knot_vector,iknots(i));
        D = D*Di;
    end
else
    if size(knot_vector,2) == 2 % surface
        % compute refinement matrix on the u direction
        D1 = eye(nurbs.number(1));
        for i = 1:length(iknots{1})
            [Di,knot_vector{1}] = bspkntins_coeffs_1(nurbs.order(1)-1,knot_vector{1},iknots{1}(i));
            D1 = D1*Di;
        end
        
        % compute refinement matrix on the v direction
        D2 = eye(nurbs.number(2));
        for i = 1:length(iknots{2})
            [Di,knot_vector{2}] = bspkntins_coeffs_1(nurbs.order(2)-1,knot_vector{2},iknots{2}(i));
            D2 = D2*Di;
        end
        D1
        D2
        % compute the bivariate coefficient matrix
        D = outer_prod(D2,D1);
    elseif size(knot_vector,2) == 3 % volume
        % compute refinement matrix on the u direction
        D1 = eye(nurbs.number(1));
        for i = 1:length(iknots{1})
            [Di,knot_vector{1}] = bspkntins_coeffs_1(nurbs.order(1)-1,knot_vector{1},iknots{1}(i));
            D1 = D1*Di;
        end

        % compute refinement matrix on the v direction
        D2 = eye(nurbs.number(2));
        for i = 1:length(iknots{2})
            [Di,knot_vector{2}] = bspkntins_coeffs_1(nurbs.order(2)-1,knot_vector{2},iknots{2}(i));
            D2 = D2*Di;
        end

        % compute refinement matrix on the w direction
        D3 = eye(nurbs.number(3));
        for i = 1:length(iknots{3})
            [Di,knot_vector{3}] = bspkntins_coeffs_1(nurbs.order(3)-1,knot_vector{3},iknots{3}(i));
            D3 = D3*Di;
        end

        % compute the trivariate coefficient matrix
        D = outer_prod(D3,outer_prod(D2,D1));
    else
        error('Invalid size of inserted knot vector');
    end
end


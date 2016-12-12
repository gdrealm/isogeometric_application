function varargout = nrbinterp(Pts,params)
%
% INTERP: Create a NURBS curve/surface/volume based on control points and order
% Calling Sequences:
%   new_nurbs = nrbinterp(nurbs, method, tt)

% INPUT:
%
%   params: structure contain the information controlling the interpolation
%   p: order of the NURBS patch
%

p = params.p;
interp_method = params.interp_method;
knots_method = params.knots_method;

p_size = size(p(:));

if p_size(1) == 1
%% NURBS curve

    % bound check
    N = size(Pts, 1);
    if N < p+1
        error('Not enough control points to build the curve');
    end
    
    if strcmp(interp_method, 'equally spaced')
        % build the knot for the sampling points; assuming that the knots used for interpolation are equally distanced; assuming weight = 1
        interp_knots = linspace(0,1,N);
    elseif strcmp(interp_method, 'chord length')
        % build the knot for the sampling points; using the chord length
        % method (eq 9.4, The NURBS book, Les Piegl)
        d = 0.0;
        for i = 1:N-1
            d = d + norm(Pts(i+1,:) - Pts(i,:));
        end
        interp_knots = zeros(1,N);
        for i = 1:N-1
            interp_knots(i+1) = interp_knots(i) + ...
                norm(Pts(i+1,:) - Pts(i,:))/d;
        end
        interp_knots(N) = 1.0; % it is one as construction; but to avoid numerical error, I set it to one explicitly
    elseif strcmp(interp_method, 'centripetal')
        % build the knot for the sampling points; using the centripetal
        % method (eq 9.6, The NURBS book, Les Piegl)
        d = 0.0;
        for i = 1:N-1
            d = d + sqrt(norm(Pts(i+1,:) - Pts(i,:)));
        end
        interp_knots = zeros(1,N);
        for i = 1:N-1
            interp_knots(i+1) = interp_knots(i) + ...
                sqrt(norm(Pts(i+1,:) - Pts(i,:)))/d;
        end
        interp_knots(N) = 1.0; % it is one as construction; but to avoid numerical error, I set it to one explicitly
    else
        error(['Unknown interpolation method ' interp_method]);
    end
    
    if strcmp(knots_method, 'equally spaced')
        % build the knot vector
        knots = [zeros(1,p) linspace(0,1,N-p+1) ones(1,p)];
    elseif strcmp(knots_method, 'averaging')
        % build the knot vector using averaging technique (eq 9.8, The
        % NURBS book, Les Piegl)
        % TODO check this
        knots = [zeros(1,p+1) zeros(1,N-p-1) ones(1,p+1)];
        for i = 2:N-p
            tmp = 0;
            for j = i:i+p-1
                tmp = tmp + interp_knots(j);
            end
            knots(i+p) = tmp/p;
        end
    else
        error(['Unknown knots contruction method ' knots_method]);
    end
    
    % build the interpolation matrix
    A = zeros(N,N);
    for i = 1:N
        span = findspan(N-1,p,interp_knots(i),knots);
        F = basisfun(span,interp_knots(i),p,knots);
        start = span - p + 1;
        A(i,start:(start+p)) = F;
    end
        
    % compute the control points
    CtrlPts = zeros(N,3);
    for i = 1:3
        b = Pts(:,i);
        CtrlPts(:,i) = A\b;
    end
    
    % create the nurbs
    nurbs = nrbmak(CtrlPts',knots);

elseif p_size(1) == 2
%% NURBS surface
% TODO
error('nrbinterp is not implemented yet for surface');

elseif p_size(1) == 3
%% NURBS volume
% TODO
error('nrbinterp is not implemented yet for volume');

else
    error('something wrong with the order');
end

if nargout == 0
    nrbplot(nurbs, params.resolution);
elseif nargout == 1
    varargout{1} = nurbs;
elseif nargout == 2
    varargout{1} = nurbs;
    varargout{2} = interp_knots;
end

function new_nurbs = nrbslice(nurbs, method, tt)
%
% SLICE: Evaluate the slice over a NURBS surface, or volume
% Calling Sequences:
%   new_nurbs = nrbslice(nurbs, method, tt)

% INPUT:
%
%   nurbs: NURBS volume or surface
%
%   method: 'u'     slice in u direction => for volume, it returns the surface
%                sliced in 'vw'; for surface, it returns the line sliced 
%                in 'v'; for line, it returns a point.
%   method: 'uv'    slice in uv direction => for volume, it returns the line
%                sliced in 'w'; for surface, it returns a point; for line,
%                it returns an error.
%   method: 'uvw'   slice in uvw direction => for volume, it returns a
%                point; for line and surface, it returns an error.
%   tt: array of slicing coordinate, it has the same length as number of
%       character of method; the values must also come in the same sequence of method

if (~isstruct(nurbs))
    error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
    error('Not a recognised NURBS representation');
end

if (iscell(nurbs.knots))
    if (size(nurbs.knots,2) == 2)
        %% NURBS structure represents a surface
        num1 = nurbs.number(1);
        num2 = nurbs.number(2);
        degree = nurbs.order-1;

        if (strcmp(method, 'u'))
            val = permute(nurbs.coefs,[1 3 2]);
            val = reshape(val, 4*num2, num1);
            val = bspeval(degree(1), val, nurbs.knots{1}, tt);
            val = reshape(val, [4 num2 1]);
            new_nurbs = nrbmak(val, nurbs.knots{2});

        elseif (strcmp(method, 'v'))
            val = reshape(nurbs.coefs, 4*num1, num2);
            val = bspeval(degree(2), val, nurbs.knots{2}, tt);
            val = reshape(val,[4 num1 1]);
            new_nurbs = nrbmak(val, nurbs.knots{1});

        elseif (strcmp(method, 'uv'))
            [p,w] = nrbeval(nurbs, tt);
            new_nurbs = [p;w];

        else
            error('Not a recognised method for slicing')
        end
            
    elseif (size(nurbs.knots,2) == 3)
        %% NURBS structure represents a volume
        num1 = nurbs.number(1);
        num2 = nurbs.number(2);
        num3 = nurbs.number(3);
        degree = nurbs.order-1;

        if (strcmp(method, 'u'))
            val = permute(nurbs.coefs, [1 3 4 2]); %[1 3 2 4](possibly wrong)
            val = reshape(val, 4*num2*num3, num1);
            val = bspeval(degree(1), val, nurbs.knots{1}, tt);
            val = reshape(val, [4 num2 num3 1]);
            new_nurbs = nrbmak(val, {nurbs.knots{2} nurbs.knots{3}});

        elseif (strcmp(method, 'v'))
            val = permute(nurbs.coefs,[1 2 4 3]); %[1 2 4 3] %[1 4 2 3]
            val = reshape(val, 4*num1*num3, num2);
            val = bspeval(degree(2), val, nurbs.knots{2}, tt);
            val = reshape(val, [4 num1 num3 1]);
            new_nurbs = nrbmak(val, {nurbs.knots{1} nurbs.knots{3}});

        elseif (strcmp(method, 'w'))
            val = reshape (nurbs.coefs, 4*num1*num2, num3);
            val = bspeval (degree(3), val, nurbs.knots{3}, tt);
            val = reshape (val, [4 num1 num2 1]);
            new_nurbs = nrbmak(val, {nurbs.knots{1} nurbs.knots{2}});

        elseif (strcmp(method, 'uv'))
            srf = nrbslice(nurbs, 'u', tt{1});
            new_nurbs = nrbslice(srf, 'u', tt{2});

        elseif (strcmp(method, 'vw'))
            srf = nrbslice(nurbs, 'v', tt{1});
            new_nurbs = nrbslice(srf, 'v', tt{2});

        elseif (strcmp(method, 'uw'))
            srf = nrbslice(nurbs, 'u', tt{1});
            new_nurbs = nrbslice(srf, 'v', tt{2});

        elseif (strcmp(method, 'uvw'))
            [p,w] = nrbeval(nurbs, tt);
            new_nurbs = [p;w];

        else
            error('Not a recognised method for slicing')
        end
        
    end
else
    %% NURBS structure represents a curve
    if (strcmp(method, 'u'))
        [p,w] = nrbeval(nurbs, tt);
        new_nurbs = [p;w];
    else
        error('Not a recognised method for slicing')
    end
end


end

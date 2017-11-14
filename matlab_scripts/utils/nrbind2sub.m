%% convert the total index of the basis function to sub-index in each direction
function y = nrbind2sub(nurbs,fun)
if ~iscell(nurbs.knots) % curve
    if ind > nurbs.number
        error('The input index vector is invalid');
    else
        y = fun;
    end
else
    if size(nurbs.knots,2) == 2 % surface
        n1 = nurbs.number(1);
        n2 = nurbs.number(2);
        if fun > n1*n2 || fun < 0
            error('The index is out of range');
        end
        i2 = mod(fun-1,n2)+1;
        i1 = (fun-i2)/n2+1;
        y = [i1 i2];
    elseif size(nurbs.knots,2) == 3 % volume
        n1 = nurbs.number(1);
        n2 = nurbs.number(2);
        n3 = nurbs.number(3);
        if fun > n1*n2*n3 || fun < 0
            error('The index is out of range');
        end
        i3 = mod(fun-1,n3)+1;
        i21 = (fun-i3)/n3+1;
        i2 = mod(i21-1,n2)+1;
        i1 = (i21-i2)/n2+1;
        y = [i1 i2 i3];
    else
        error('Invalid knot vector of NURBS');
    end
end


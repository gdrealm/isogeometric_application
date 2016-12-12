%% convert the sub-index to total index 
function y = nrbsub2ind(nurbs,ind)
if ~iscell(nurbs.knots) % curve
    if ind > nurbs.number
        error('The input index vector is invalid');
    else
        y = ind;
    end
else
    if size(nurbs.knots,2) == 2 % surface
        if length(ind) ~= 2
            error('The index must be of dimension 2 for surface');
        end
        i1 = ind(1);
        i2 = ind(2);
        if i1 > nurbs.number(1) || i1 < 0
            error('Sub-index 1 is out of range');
        end
        if i2 > nurbs.number(2) || i2 < 0
            error('Sub-index 2 is out of range');
        end
        y = (i1-1) * nurbs.number(2) + i2;
    elseif size(nurbs.knots,2) == 3 % volume
        if(length(ind)~=3)
            error('The index must be of dimension 3 for volume');
        end
        i1 = ind(1);
        i2 = ind(2);
        i3 = ind(3);
        if i1 > nurbs.number(1) || i1 < 0
            error('Sub-index 1 is out of range');
        end
        if i2 > nurbs.number(2) || i2 < 0
            error('Sub-index 2 is out of range');
        end
        if i3 > nurbs.number(3) || i3 < 0
            error('Sub-index 3 is out of range');
        end
        y = ((i1-1) * nurbs.number(2) + i2 - 1) * nurbs.number(3) + i3;
    else
        error('Invalid knot vector of NURBS');
    end
end


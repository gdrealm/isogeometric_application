function y = compare_coefs(nurbs1,nurbs2,tol)
    s1 = length(nurbs1.coefs(:));
    s2 = length(nurbs2.coefs(:));
    if s1 ~= s2
        y = 1;
        return;
    end
    if norm(nurbs1.coefs(:) - nurbs2.coefs(:)) < tol
        y = 0;
    else
        y = 1;
    end
end


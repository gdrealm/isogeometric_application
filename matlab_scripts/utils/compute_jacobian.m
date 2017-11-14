function jac = compute_jacobian(nurbs,tt)
dnurbs = nrbderiv(nurbs);
[pnt, jac] = nrbdeval(nurbs,dnurbs,tt);
end


function y = h_refine(nurbs, desired_num_knots, desired_continuity)

[rknots, zeta, nknots] = kntrefine(nurbs.knots, desired_num_knots, nurbs.order - 1, desired_continuity);
new_nurbs = nrbkntins(nurbs, nknots);

y = new_nurbs;

%example
%new_nurbs = h_refine(nurbs, [2 2], [0 0])

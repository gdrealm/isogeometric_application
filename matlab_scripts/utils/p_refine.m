function y = p_refine(nurbs, desired_order)

degelev = max(desired_order - (nurbs.order - 1), 0);
new_nurbs = nrbdegelev(nurbs, degelev);

y = new_nurbs;

%example:
%new_nurbs = p_refine(nurbs, [3 3])

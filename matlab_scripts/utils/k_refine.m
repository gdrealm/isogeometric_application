function y = k_refine(nurbs, desired_order)

degelev = max(desired_order - (nurbs.order - 1), 0);
new_nurbs = nrbdegelev(nurbs, degelev);

y = new_nurbs;

% % new_nurbs.knots{1}
% % new_nurbs.knots{2}
% % new_nurbs.order - 2
% 
% 
% %     knots:      initial knot vector.
% %     n_sub:      number of new knots to be added in each interval.
% %     degree:     polynomial degree of the refined knot vector
% %     regularity: maximum global regularity 
% [rknots, zeta, nknots] = kntrefine(new_nurbs.knots, ...
%     ones(1, size(nurbs.order, 2)), ...
%     new_nurbs.order - 1, new_nurbs.order - 2) ;
% 
% % nknots{1}
% % nknots{2}
% % nknots{3}
% new_nurbs = nrbkntins(new_nurbs, nknots);
% 
% y = new_nurbs;

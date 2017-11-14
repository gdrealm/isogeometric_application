%% Plot the knot vector in bar graph
function plot_knot_vector(U)
X = unique(U);
Y = histc(U,X);
bar(X,Y,'yellow');
% bar_h = bar(X,Y);
% bar_child=get(bar_h,'Children');
% set(bar_child,'CData',sin(X));
end


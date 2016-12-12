%% subroutine to write data file to plot by gnuplot
function export_gnuplot(fn,x,y)

fid = fopen(fn, 'wt');
n = length(x);
if n~=length(y)
    error('Invalid input');
end

for i = 1:n
    fprintf(fid,'%f %f\n',x(i),y(i));
end

fclose(fid);


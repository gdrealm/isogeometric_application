%% plot a rectangle patch
% X, Y: array of x,y-coordinates
% z: the z-coordinate to plot the patch
function plot_rectangle_patch(X,Y,z,color)

verts = zeros(length(X) * length(Y),3);
for j = 1:length(Y)
    for i = 1:length(X)
        row = (j-1)*length(X) + i;
        verts(row,1) = X(i);
        verts(row,2) = Y(j);
        verts(row,3) = z;
    end
end

faces = zeros((length(X)-1)*(length(Y)-1),4);
for j = 1:length(Y)-1
    for i = 1:length(X)-1
        row = (j-1)*(length(X)-1) + i;
        start = (j-1)*length(X) + i;
        faces(row,:) = [start start+1 start+1+length(X) start+length(X)];
    end
end

p = patch('Faces',faces,'Vertices',verts,'FaceColor',color);

end


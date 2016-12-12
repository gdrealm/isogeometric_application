%% plot a box patch
% X, Y, Z: array of x,y-coordinates
function plot_box_patch(X,Y,Z,color,alpha)

verts = zeros(length(X) * length(Y) * length(Z),3);
for k = 1:length(Z)
    for j = 1:length(Y)
        for i = 1:length(X)
            row = getindex(i,j,k,length(X),length(Y),length(Z));
            verts(row,1) = X(i);
            verts(row,2) = Y(j);
            verts(row,3) = Z(k);
        end
    end
end

for i = 1:length(X)
    faces = zeros((length(Y)-1)*(length(Z)-1),4);
    row = 1;
    for k = 1:length(Z)-1
        for j = 1:length(Y)-1
            n1 = getindex(i,j,k,length(X),length(Y),length(Z));
            n2 = getindex(i,j+1,k,length(X),length(Y),length(Z));
            n3 = getindex(i,j+1,k+1,length(X),length(Y),length(Z));
            n4 = getindex(i,j,k+1,length(X),length(Y),length(Z));
            faces(row,:) = [n1 n2 n3 n4];
            row = row + 1;
        end
    end
    patch('Faces',faces,'Vertices',verts,'FaceColor',color,'FaceAlpha',alpha);
end

for j = 1:length(Y)
    faces = zeros((length(X)-1)*(length(Z)-1),4);
    row = 1;
    for i = 1:length(X)-1
        for k = 1:length(Z)-1
            n1 = getindex(i,j,k,length(X),length(Y),length(Z));
            n2 = getindex(i+1,j,k,length(X),length(Y),length(Z));
            n3 = getindex(i+1,j,k+1,length(X),length(Y),length(Z));
            n4 = getindex(i,j,k+1,length(X),length(Y),length(Z));
            faces(row,:) = [n1 n2 n3 n4];
            row = row + 1;
        end
    end
    patch('Faces',faces,'Vertices',verts,'FaceColor',color,'FaceAlpha',alpha);
end

for k = 1:length(Z)
    faces = zeros((length(X)-1)*(length(Y)-1),4);
    row = 1;
    for j = 1:length(Y)-1
        for i = 1:length(X)-1
            n1 = getindex(i,j,k,length(X),length(Y),length(Z));
            n2 = getindex(i+1,j,k,length(X),length(Y),length(Z));
            n3 = getindex(i+1,j+1,k,length(X),length(Y),length(Z));
            n4 = getindex(i,j+1,k,length(X),length(Y),length(Z));
            faces(row,:) = [n1 n2 n3 n4];
            row = row + 1;
        end
    end
    patch('Faces',faces,'Vertices',verts,'FaceColor',color,'FaceAlpha',alpha);
end

end

function y = getindex(i,j,k,lenx,leny,lenz)
    y = ((k-1)*leny + j-1)*lenx + i;
end


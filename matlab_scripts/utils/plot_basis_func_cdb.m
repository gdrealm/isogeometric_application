%% plot the NURBS basis function using Cox-de-Boor formula
function plot_basis_func_cdb(funs,Id,Xi,Eta,params)
num_points1 = params.num_points1;
num_points2 = params.num_points2;
p1 = params.p1;
p2 = params.p2;

for k = 1:length(funs)
    for l = 1:length(Id)
        if funs(k) == Id(l)
            Xi_local = Xi{l};
            Eta_local = Eta{l};
        end
    end
    x_xi = linspace(min(Xi_local),0.999999*max(Xi_local),num_points1);
    x_eta = linspace(min(Eta_local),0.999999*max(Eta_local),num_points2);
    X = zeros(num_points1,num_points2);
    Y = zeros(num_points1,num_points2);
    Z = zeros(num_points1,num_points2);
    for i = 1:num_points1
        for j = 1:num_points2
            X(i,j) = x_xi(i);
            Y(i,j) = x_eta(j);
            Z(i,j) = CoxDeBoor(x_xi(i),1,p1,Xi_local) * CoxDeBoor(x_eta(j),1,p2,Eta_local);
        end
    end
    surf(X,Y,Z,'FaceColor','yellow');
end

% params.num_points1 = 10;
% params.num_points2 = 10;
% params.p1 = 1;
% params.p2 = 2;
% funs = [955];
% plot_basis_func_cdb(funs,Id,Xi,Eta,params);

%% plot the hierarchical NURBS geometry using Cox-de-Boor formula
function A = check_linear_independence(Xi,Eta,W,Id,params)
min_xi = params.min_xi;
max_xi = params.max_xi;
min_eta = params.min_eta;
max_eta = params.max_eta;
num_points1 = params.num_points1;
num_points2 = params.num_points2;
x_xi = linspace(min_xi,max_xi,num_points1);
x_eta = linspace(min_eta,max_eta,num_points2);
p1 = params.p1;
p2 = params.p2;
A = zeros(num_points1*num_points2,length(W));
for i = 1:num_points1
    for j = 1:num_points2
        row = (i-1)*num_points2 + j;
        Denom = 0.0;
        for k = 1:length(W)
            N = W(k) * CoxDeBoor(x_xi(i),1,p1,Xi{k}) * CoxDeBoor(x_eta(j),1,p2,Eta{k});
            A(row,k) = N;
            Denom = Denom + N;
        end
        for k = 1:length(W)
            A(row,k) = A(row,k) / Denom;
        end
    end
end
disp(['assembly of matrix A completed'])
size(A)
B = null(A)
r = length(W) - size(B,2)
if(r ~= length(W))
    disp('The basis functions is not linear independent with the provided sampling points');
    dependent_bfs = [];
    for i = 1:length(W)
        if(abs(B(i,1)) > 1.0e-10)
            dependent_bfs = [dependent_bfs Id(i)];
        end
    end
    dependent_bfs
end

% params.num_points1 = 51;
% params.num_points2 = 51;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.p1 = 2;
% params.p2 = 2;
% A = check_linear_independence(Xi,Eta,W,Id,params);



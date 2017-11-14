%% plot the NURBS basis function using Cox-de-Boor formula
function plot_basis_func_cdb(funs,S,C,N,params)
num_points1 = params.num_points1;
num_points2 = params.num_points2;
p1 = params.p1;
p2 = params.p2;

for f = 1:length(funs)
    % compute and plot the basis function from Bezier decomposition
    basis_id = funs(f);
    cnt = 1;
    for i = 1:length(N)
         for j = 1:length(N{i})
             if N{i}(j) == basis_id
                 basis_in_elem(cnt,:)=[i j];
                 cnt = cnt+1;
             end
         end
    end
    for n = 1:size(basis_in_elem,1)
         id    = basis_in_elem(n,1);
         a_xi  = S{id}(1,1);
         b_xi  = S{id}(1,2);
         a_eta = S{id}(2,1);
         b_eta = S{id}(2,2);
         x_xi  = linspace(a_xi,b_xi,num_points1);
         x_eta = linspace(a_eta,b_eta,num_points2);
         X = zeros(num_points1,num_points2);
         Y = zeros(num_points1,num_points2);
         Z = zeros(num_points1,num_points2);
         t = basis_in_elem(n,2);
         Crow = C{id}(t,:);
         for i = 1:num_points1
             for j = 1:num_points2
                 X(i,j) = x_xi(i);
                 Y(i,j) = x_eta(j);
                 for k = 0:p1
                     y_xi = bezier((x_xi(i)-a_xi)/(b_xi-a_xi),k,p1);
                     for l = 0:p2
                         y_eta  = bezier((x_eta(j)-a_eta)/(b_eta-a_eta),l,p2);
                         Z(i,j) = Z(i,j) + y_xi * y_eta * Crow(k*(p2+1)+l+1);
                     end
                end
             end
         end
         surf(X,Y,Z);
    end
end

% params.num_points1 = 10;
% params.num_points2 = 10;
% params.p1 = 1;
% params.p2 = 2;
% funs = [950 955];
% plot_basis_func_bezier(funs,S,C,N,params);


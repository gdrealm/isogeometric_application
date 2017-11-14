%% plot the NURBS basis function using Cox-de-Boor formula
% funs: list of bf id to plot
% z: z-coordinate to put the plot
% color: plot color
function plot_basis_func_support_cells_at(funs,Id,Xi,Eta,z,color)
for k = 1:length(funs)
    for l = 1:length(Id)
        if funs(k) == Id(l)
            Xi_local = Xi{l};
            Eta_local = Eta{l};
            plot_rectangle_patch(Xi_local,Eta_local,z,color);
        end
    end
end


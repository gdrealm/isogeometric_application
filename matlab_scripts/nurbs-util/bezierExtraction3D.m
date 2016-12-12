function [C, Cxi, Ceta, Czeta] = bezierExtraction3D(uknot, vknot, wknot, p, q, r)
%
% Bezier extraction operators for a 3D NURBS solid.
%
% HG Bui
% Based on the code bezierExtraction2D of
% VP Nguyen
% Cardiff University, UK

% Bezier extraction operators for xi and eta, zeta
% nb1: number of elements along xi direction
% nb2: number of elements along eta direction
% nb3: number of elements along zeta direction

[Cxi, nb1]   = bezierExtraction(uknot, p);
[Ceta, nb2]  = bezierExtraction(vknot, q);
[Czeta, nb3] = bezierExtraction(wknot, r);

% For Bsplines/NURBS, the element Bezier extraction
% operator is square.

size1 = size(Cxi(:, :, 1), 1);
size2 = size(Ceta(:, :, 1), 1);
size3 = size(Czeta(:, :, 1), 1);

% Bezier extraction operators for the whole mesh
% as the trivariate tensor product of Cxi, Ceta and Czeta
C_eta_xi = zeros(size1 * size2, size1 * size2);
C = zeros(size1 * size2 * size3, size1 * size2 * size3, nb1 * nb2 * nb3);

for eta = 1:nb2
    for xi = 1:nb1
        %firstly compute C^j_\eta \outer_prod C^k_\xi
        for row = 1:size2
            ird = (row - 1) * size1 + 1;
            jrd =  row * size1;
            for col = 1:size2
                icd = (col - 1) * size1 + 1;
                jcd =  col * size1;
                C_eta_xi(ird:jrd, icd:jcd) = Ceta(row, col, eta) * Cxi(:, :, xi);
            end
        end
        %secondly compute C^i_\zeta \outer_prod C_\eta\xi
        for zeta = 1:nb3
            e = ((zeta - 1) * nb2 + eta - 1) * nb1 + xi;
            for row = 1:size3
                ird = (row - 1) * size1 * size2 + 1;
                jrd =  row * size1 * size2;
                for col = 1:size3
                    icd = (col - 1) * size1 * size2 + 1;
                    jcd =  col * size1 * size2;
                    C(ird:jrd, icd:jcd, e) = Czeta(row, col, zeta) * C_eta_xi;
                end
            end
        end
    end
end


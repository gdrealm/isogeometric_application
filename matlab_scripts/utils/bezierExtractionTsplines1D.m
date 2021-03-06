% Bezier extraction operator computation algorithm for B-spline basis function with local knot vector in 1D
% Reference: Scott et al, Isogeometric finite element data structures based on Bezier extraction of T-splines
function [C,nb,Ubar] = bezierExtractionTsplines1D(Xi,U,spans,p)
% Input:
%   Xi      local knot vector
%   U       interior knots to be inserted into Xi (need to be non-repeated)
%   spans   knot spans in Xi where new knots will be inserted (need to be the same size as U)
%   p       curve degree
%
% Output:
%   C       matrix of size nb x (p+1) contain the rows of bezier extraction of the T-splines basis function at each knot span of the filled extended knot vector
%   nb      number of spans of the filled extended knot vector
%   Ubar    filled extended knot vector

% do some check
if length(Xi) ~= p+2
    error('local knot vector must be of length p + 2')
end

if length(U) ~= length(spans)
    error('spans array must be the same length as the number of interior knots')
end

[Ubar,nt] = compute_extended_knot_vector(Xi,p);
a = p + 1;
b = a + 1;
nb = 1;
C(1,:) = zeros(1,p+1);
C(1,nt+1) = 1;
m = length(U);
mbar = p + 2 + nt + m;
ki = 1;
si = 1;
total_add = 0;
while b < mbar
    % count multiplicity of knots at location b
    add = 0;
    if si <= m && spans(si) == ki
        mult = 0;
        add = 1;
        % add the new knot to the knot vector
%         Ubar(b+1:mbar+p-m+si) = Ubar(b:mbar+p-m+si-1);
        Ubar(b+1:length(Ubar)+1) = Ubar(b:length(Ubar));
        Ubar(b) = U(si);
        si = si + 1;
    else
        ki = ki + 1;
        i = b;
        b
        mbar
        Ubar
%        while b <= mbar && Ubar(b+1) == Ubar(b)
        while b < mbar
            if Ubar(b+1) == Ubar(b)
                b = b + 1;
                if b >= mbar
                    break
                end
            else
                break
            end
        end
        mult = b - i + 1;
    end
    total_add = total_add + add; % count the total number of additional knots

%     b
%     i
%     mult
    
    if mult <= p
        C(nb+1,:) = zeros(1,p+1); % initialize the next extraction operator row
        loc = nt+1-nb+total_add; % identify the next location to be 1
        if loc >= 1 && loc <= p+1
            C(nb+1,loc) = 1.0;
        end
    
        numer = Ubar(b) - Ubar(a);
        for j = p:-1:mult+1
            alphas(j - mult) = numer / (Ubar(a + j + add) - Ubar(a));
        end
        r = p - mult;
        % update the matrix coefficients for r new knots
        for j = 1:r
            save = r - j + 1;
            s = mult + j;
            for k = (p+1):-1:(s+1)
                alpha = alphas(k - s);
                C(nb,k) = alpha * C(nb,k) + (1 - alpha) * C(nb,k-1);
            end
            if b <= mbar
                % update overlapping coefficients of the next operator row
                C(nb+1,save) = C(nb,p+1);
            end
        end
        nb = nb + 1; % finished with the current operator
        if b <= mbar
            % update indices for the next operator
            a = b;
            b = b + 1;
        end
%        disp(['-----end of branch mult < p']);
%     elseif mult == p
%         if b <= mbar
%             nb = nb + 1;
%             a = b;
%             b = b + 1;
%         end
%        disp(['-----end of branch mult == p']);
    end
end

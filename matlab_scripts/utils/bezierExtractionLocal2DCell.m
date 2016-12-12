function [C,nt_xi,nt_eta,Ubar_xi,Ubar_eta] = bezierExtractionLocal2DCell(Xi,Eta,Uxi,Ueta,p,q)

[Cxi,nt_xi,Ubar_xi]    = bezierExtractionLocal1DCell(Xi,Uxi,p);
[Ceta,nt_eta,Ubar_eta] = bezierExtractionLocal1DCell(Eta,Ueta,q);

C = zeros((p+1)*(q+1));

for i = 1:p+1
    for j = 1:q+1
        C((i-1)*(q+1)+j) = Cxi(i)*Ceta(j);
    end
end


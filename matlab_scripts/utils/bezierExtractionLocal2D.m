function [C,nb_xi,nb_eta,Ubar_xi,Ubar_eta] = bezierExtractionLocal2D(Xi,Eta,Uxi,Ueta,p,q)

[Cxi,nb_xi,Ubar_xi]    = bezierExtractionLocal1D(Xi,Uxi,p)
[Ceta,nb_eta,Ubar_eta] = bezierExtractionLocal1D(Eta,Ueta,q)

C = zeros(nb_xi*nb_eta,(p+1)*(q+1));

for i = 1:size(Cxi,1)
    for j = 1:size(Cxi,2)
        C((i-1)*size(Ceta,1)+1 : i*size(Ceta,1),...
            (j-1)*size(Ceta,2)+1 : j*size(Ceta,2)) = Cxi(i,j)*Ceta;
    end
end

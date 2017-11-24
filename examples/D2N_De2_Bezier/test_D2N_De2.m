addpath('/home/hbui/workspace2/isogeometric/nurbs');

knts = [0 0 0 1 1 1];

xi = 0.112702;
eta = 0.112702;

p = 2;
n = 6-p-1;

sxi = findspan(n-1,p,xi,knts);
NDxi = basisfunder(s,p,xi,knts,2);
Nxi = reshape(NDxi(1,1,:),1,3);
Nxi1 = reshape(NDxi(1,2,:),1,3);
Nxi2 = reshape(NDxi(1,3,:),1,3);

seta = findspan(n-1,p,eta,knts);
NDeta = basisfunder(s,p,eta,knts,2);
Neta = reshape(NDeta(1,1,:),1,3);
Neta1 = reshape(NDeta(1,2,:),1,3);
Neta2 = reshape(NDeta(1,3,:),1,3);

dNdxi = zeros(1,9);
dNdeta = zeros(1,9);
d2Ndxi2 = zeros(1,9);
d2Ndeta2 = zeros(1,9);
d2Ndxideta = zeros(1,9);
for i = 1:3
    for j = 1:3
        k = 3*(j-1)+i;
        dNdxi(k) = Nxi1(i)*Neta(j);
        dNdeta(k) = Nxi(i)*Neta1(j);
        d2Ndxi2(k) = Nxi2(i)*Neta(j);
        d2Ndxideta(k) = Nxi1(i)*Neta1(j);
        d2Ndeta2(k) = Nxi(i)*Neta2(j);
    end
end

Nxi'*Neta

dNdxi
dNdeta

d2Ndxi2
d2Ndeta2
d2Ndxideta


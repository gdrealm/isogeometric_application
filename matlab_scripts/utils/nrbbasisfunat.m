%% evaluate a specific B-spline basis function based on index at the specific points
function Bf = nrbbasisfunat(points,fun,nurbs,type)
% Input:
%   u       evaluation points
%   fun     index of the basis function
%   nurbs   nurbs data structure
% Output:
%   Bf      basis function value of index fun at evaluation points

if ~iscell(nurbs.knots) % curve
    if max(fun) > nurbs.number || max(fun) > nurbs.number || min(fun) < 1 || max(fun) < 1
        error('index of basis function is out of size');
    end

    % compute the basis function and index array of nonzero functions
    [B,id] = nrbbasisfun(points,nurbs);

    %rearrange to the correct sequence
    numbasis = nurbs.number;
    Bn = zeros(size(B,1),numbasis);
    for i = 1:size(B,1)
        for j = 1:size(B,2)
            Bn(i,id(i,j)+1) = B(i,j);
        end
    end

    Bf = Bn(:,fun);
else
    if size(nurbs.knots,2) == 2 % surface
        if ~iscell(fun)
            % if fun is not given as cell, then it will be parsed to have 
            [fun1,fun2] = ind2sub([nurbs.number(1),nurbs.number(2)],fun);
        else
            fun1 = fun{1};
            fun2 = fun{2};
%             fun = sub2ind([nurbs.number(1),nurbs.number(2)],fun1,fun2);
        end

        if ~iscell(points)
            error('the points must be given by cell')
        else
            u = points{1};
            v = points{2};
        end

        % compute the basis function and index array of nonzero functions
        [B,id] = nrbbasisfun({u,v},nurbs);

        %rearrange to the correct sequence
        Bn = zeros(length(u),length(v),nurbs.number(1)*nurbs.number(2));
        for i = 1:size(B,1) %length(u)*length(v)
            for j = 1:size(B,2) %(p+1)*(q+1)
                [row,col] = ind2sub([length(u),length(v)],i);
                Bn(row,col,id(i,j)) = B(i,j);
            end
        end

        %rearrange to the correct sequence
        if strcmp(type,'sub')
            Bf = zeros(length(u),length(v),length(fun1),length(fun2));
            for i = 1:length(fun1)
                for j = 1:length(fun2)
                    k = sub2ind([nurbs.number(1),nurbs.number(2)],fun1(i),fun2(j));
                    Bf(:,:,i,j) = Bn(:,:,k);
                end
            end
        elseif strcmp(type,'ind')
            Bf = Bn;
        else
            error('Invalid type of output function');
        end
    elseif size(nurbs.knots,2) == 3 % volume
        % TODO
    else
        error('Invalid knot vector of NURBS');
    end
end


%% example usage
% coefs = [
%    0.00  1.00  2.00  3.00  4.00  5.00  6.00
%    0.00  0.20  1.00  0.30  0.40  0.20  0.10
%    0.00  0.00  0.00  0.00  0.00  0.00  0.00
%    1.00  1.00  1.00  1.00  1.00  1.00  1.00];
% 
% knots = [0 0 0 0 1 2 3 4 4 4 4];
% 
% nurbs = nrbmak(coefs,knots);
% 
% fun = [1 2 3 4 5 6 7];
% n = 100;
% u = linspace(min(nurbs.knots),max(nurbs.knots),n);
% B = nrbbasisfunat(u,fun,nurbs);
% plot(u,B);


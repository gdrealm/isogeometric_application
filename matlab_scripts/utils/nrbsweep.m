function surf = nrbsweep(curve,traj,v)

% 
% nrbsweep: Construct a NURBS surface by sweeping a NURBS curve along
% a NURBS curve, or construct a NURBS volume by sweeping a NURBS surface
% along a NURBS curve; using the skinning technique (Algorithm A10.2 The NURBS book, Les Pielg)
%
% Calling Sequence:
% 
%   srf = nrbsweep(crv,traj,v)
% 
% INPUT:
% 
%   crv		: NURBS curve or surface to revolve, see nrbmak.
% 
%   traj	: trajectory of the sweep. Need to be a NURBS curve.
%
%   v       : sampling knot on trajectory to create the skinning
% 
% OUTPUT:
%
%   srf		: constructed surface or volume
% 
% Description:
% 
%   Construct a NURBS surface by sweeping the profile NURBS curve along
%   an axis defined by a NURBS curve.
% 
% Examples:
%       TODO
%
% NOTE:
%
%   The algorithm:
%       TODO
%

if (nargin < 2)
  error('Not enough arguments to construct sweeped surface');
end

if (iscell(curve.knots) && numel(curve.knots) == 3)
  error('The function nrbsweep is not yet ready to create volumes') 
end

if (iscell (curve.knots))
    %% Construct the sweeped volume
    % compute the frame at each section
    n = size(v(:),1);
    sections = cell(n);
    for i = 1:n
        % curve section point, first and second derivative
        [ders,ders2] = nrbderiv(traj);
        [o,d,d2] = nrbdeval(traj,ders,ders2,v(i));
        x = d/norm(d);
        
        % this method of computing B can cause twist of the swept volume
%         aux = cross(d,d2);
%         B = aux/norm(aux);
        
        if i == 1
            aux = cross([1 0 0]',x);
            B = aux/norm(aux);
        else
            aux = B - dot(B,x)*x;
            B = aux/norm(aux);
        end

        z = B/norm(B);
        y = cross(z,x);
        A = [x' 0;y' 0; z' 0;o' 1]';
        if i > 1
            new_curve = nrbtform(curve,A/A1);
            sections{i} = new_curve;
        else
            sections{i} = curve;
            A1 = A;
        end
    end
    
%     for i=1:n
%         plot_ctrl_points_2d(sections{i},'off');
%     end
    
    % for each point in the surface create the interpolation curve
    interp_params.p = traj.order - 1;
    interp_params.interp_method = 'centripetal';
    interp_params.knots_method = 'averaging';
    Pts = zeros(n,4);
    interp_curves = cell(size(curve.coefs,2),size(curve.coefs,3));
    for i = 1:size(curve.coefs,2)
        for j = 1:size(curve.coefs,3)
            for k = 1:n
                Pts(k,:) = sections{k}.coefs(:,i,j)';
            end
            interp_curves{i,j} = nrbinterp(Pts,interp_params);
%             interp_curve
%             plot_ctrl_points_1d(interp_curve,'off');
%             interp_curves{i,j}.knots
        end
    end
    surf_knots{1} = curve.knots{1};
    surf_knots{2} = curve.knots{2};
    surf_knots{3} = interp_curves{1,1}.knots;
    
    surf_coefs = zeros(4,size(curve.coefs,2),size(curve.coefs,3),size(interp_curves{1,1}.coefs,2));
    
    for i = 1:size(curve.coefs,2)
        for j = 1:size(curve.coefs,3)
            for k = 1:size(interp_curves{1,1}.coefs,2)
                surf_coefs(:,i,j,k) = interp_curves{i,j}.coefs(:,k);
            end
        end
    end
    
    surf = nrbmak(surf_coefs, surf_knots);
else
    %% Construct the sweeped surface

end

end


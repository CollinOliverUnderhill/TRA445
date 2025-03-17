function dcdt = diffusion_sphere_withGhost_var(c, R, D, N)
    % diffusion_sphere_withGhost_var:
    %   radial diffusion in spherical coords,
    %   uses ghost node approach at r=0.

    dcdt = zeros(size(c));
    dr   = R / N;
    r    = linspace(0, R, N+1)';

    % r=0 boundary (ghost node)
    i = 1;
    c_ghost = c(i+1);  % approximate symmetry => c_ghost = c(2)
    d2dr2_0 = (c_ghost - 2*c(i) + c_ghost)/(dr^2);
    dcdt(i) = D * d2dr2_0;

    % interior nodes
    for i = 2:N
        d2dr2 = (c(i+1) - 2*c(i) + c(i-1))/dr^2;
        dcdr  = (c(i+1) - c(i-1)) / (2*dr);
        dcdt(i) = D * ( d2dr2 + (2/r(i))*dcdr );
    end
    % outer boundary (i=N+1) flux correction is done externally
end



% function dcdt = diffusion_sphere_withGhost_var(c, R, D, N)
%     % c: vector of length N+1
%     % R: particle radius
%     % D_ref: diffusion coefficient
%     % N: # of sub-intervals in radial discretization
% 
%     dr = R / N;
%     r  = linspace(0, R, N+1)';
% 
%     dcdt = zeros(size(c));
% 
%     % Interior nodes
%     for i = 2:N
%         dcdt(i) = D * ( ...
%             (c(i+1) - 2*c(i) + c(i-1)) / dr^2 + ...
%             (2/r(i)) * (c(i+1) - c(i-1)) / (2*dr) );
%     end
%     dcdt(N+1)=D*((2*c(N)-2*c(N+1))/dr^2 );
%     % Boundary at r=0 (due to symmetry)
%     dcdt(1) = 2 * D * ...
%     (c(2) - c(1)) / dr^2;  % or a specialized second-order approach
% end
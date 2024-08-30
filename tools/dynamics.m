function [xp,dfdx,dfdw] = dynamics(x,Deltax,gyr,dt)
% DYNAMICS - Propagate state through dynamic model and obtain matrices for
% covariance update
%
% Syntax:
%   [xp,dfdx,dfdw] = dynamics(x,Deltax,gyr,dt)
%
% In:
%   x           - State
%   Deltax      - Delta position measurements
%   gyr         - Angular velocity measurements
%   dt          - Time between samples
%
% Out:
%   xp          - Predicted state
%   dfdx        - Jacobian wrt state
%   dfdw        - Jacobian wrt noise
%
% Description:
%   Propagate state through dynamics to get prediction and Jacobian
%   matrices. See [1] for details.
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Update state
% Copy state to make sure that all landmark states 
% as well as the bias state are simply propagated
xp = x; 

% Update position and heading states
psi = x(3);
xp(1:2) = x(1:2) + [ cos(psi) sin(psi); 
                   -sin(psi) cos(psi)]*Deltax;
xp(3) = x(3) + (gyr-x(4))*dt;

%% Compute Jacobians
% Jacobian w.r.t. state itself
dfdx = eye(numel(x));
dfdx(3,4) = -dt; 
dfdx(1:2,3) = [-sin(psi)  cos(psi);
               -cos(psi) -sin(psi)]*Deltax;

% Jacobian w.r.t. state dynamic noise
dfdw = eye(numel(x),3); 
dfdw(1:2,1:2) = [ cos(psi) sin(psi);  
                   -sin(psi) cos(psi)];
dfdw(3,3) = dt;

end
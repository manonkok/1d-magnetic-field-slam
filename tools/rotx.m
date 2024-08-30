function R=rotx(theta)

% ROTX - Rotation matrix for a rotation around the x-axis by an angle (radians)
%
% Syntax:
%   R=rotx(theta)
%
% In:
%   theta           - Angle (radians)
%
% Out:
%   R               - Rotation matrix
%
% Description:
%   Rotation matrix for a rotation around the x-axis by an angle (radians)
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Compute rotation matrix
R = [1 0 0 ; ...
    0 cos(theta) sin(theta) ; ...
    0, -sin(theta), cos(theta)];

end

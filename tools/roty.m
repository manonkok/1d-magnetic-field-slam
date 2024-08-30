function R=roty(theta)

% ROTY - Rotation matrix for a rotation around the y-axis by an angle (radians)
%
% Syntax:
%   R=roty(theta)
%
% In:
%   theta           - Angle (radians)
%
% Out:
%   R               - Rotation matrix
%
% Description:
%   Rotation matrix for a rotation around the y-axis by an angle (radians)
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Compute rotation matrix
R = [cos(theta), 0, -sin(theta) ; ...
        0, 1, 0 ; ...
        sin(theta), 0, cos(theta)];

end


function e = quat2eul(q,varargin)
%% QUAT2EUL - Converts quaternions into Euler angles
%
% Syntax:
%   e = quat2eul(q,seq)
%
% In:
%   q   - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%         both 1x4 and 4x1 formats are accepted.
%   seq - Sequence order ('ZYX' supported here)
%
% Out:
%   e   - Euler angles corresponding to the quaternions, in degrees. The 
%         angles are in the order of yaw, pitch, and roll.
%
% Description:
%   This function converts one or more quaternions into their corresponding 
%   Euler angles. Euler angles are a method for representing the orientation 
%   of a rigid body with a sequence of three rotations around different 
%   axes. The output is in degrees and represents yaw, pitch, and roll angles.
%
% Examples:
%   e = quat2eul([0.7071, 0, 0.7071, 0]);
%

seq = varargin{:};

if any(size(q) == 1)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
else
    q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);
end

% The parsed sequence order
switch seq
    case 'ZYX'
        as = -2*(q1.*q3-q0.*q2);
        as(as > 1) = 1;
        as(as < -1) = -1;
        
        e = [ atan2( 2*(q1.*q2+q0.*q3), q0.^2 + q1.^2 - q2.^2 - q3.^2 ), ...
            asin( as ), ...
            atan2( 2*(q2.*q3+q0.*q1), q0.^2 - q1.^2 - q2.^2 + q3.^2 )];

    otherwise
        error('Not implemented for this sequence order.')

end


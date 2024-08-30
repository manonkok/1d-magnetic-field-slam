function [dp,omega,m_b,t,pos_gt,pos_odo] = prepareData(filename,driftNoiseParams,visualiseOutput)

% PREPAREDATA - Preprocess ARKit data to use in
% one-dimensional magnetic field SLAM algorithm
%
% Syntax:
%   [dp,omega,t,pos_gt,pos_odo,psi_gt,m_b,m_w,eul] = 
%       prepareDataForFilter(filename,driftNoiseParams,visualiseOutput)
%
% In:
%   filename            - Filename to load data
%   driftNoiseParams    - Struct indicating what bias and what noise
%                           variances to use for generation of odometry data
%   visualiseOutput     - Flag if plotting ground truth and odometry paths
%
% Out:
%   dp                  - Generated delta position measurements
%   omega               - Generated angular velocity measurements
%   m_b                 - Magnetic field measurements in body frame
%   t                   - Sampling times
%   pos_gt              - Position derived from ARKit
%   pos_odo             - Position from odometry only
%   psi_gt              - Heading angle derived from ARKit
%
% Description:
%   Obtain odometry and magnetic field measurements in body frame from
%   ARKit output and add noise and bias
%   See [1] for details.
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Load data set
data = load(filename);

%% Extract position and orientation data from files
% Extract / compute pos and quat data from files
quat = data.quat;
pos_gt = [zeros(1,3) ; cumsum(data.dPos)];
    
% Cut off some additional data in the beginning
if strcmp(filename,'data/library.mat')
    indCutOff = 150;
else
    indCutOff = 1;
end

quat = quat(indCutOff:end,:);
pos_gt = pos_gt(indCutOff:end,:);
data.t = data.t(indCutOff:end,:);

% Shift data to start from origin and time to start from 0
pos_gt = pos_gt - pos_gt(1,:); % shift data to again start at (0,0)
data.t = data.t - data.t(1);
t = data.t;

% Compute Euler angles ZYX to be able to compute magnetic field in b-frame 
% and to obtain "ground truth" heading angle
eul = quat2eul(quat,'ZYX');
eul = unwrap(eul);

%% Create odometry
% Drift settings
bias = driftNoiseParams.bias; 
sh2 = driftNoiseParams.sh2; 
sp2 = driftNoiseParams.sp2; 

% Compute noise- and bias-free delta pos and delta heading
dpsi = diff(-eul(:,1));
omega = diff(-eul(:,1))./diff(t);
% Rotate delta positions from navigation to body frame
dp = diff(pos_gt(:,1:2),1);
for k=1:size(dp,1)
    tmp = sum(dpsi(1:k-1));
    dp(k,:) = ([cos(tmp) sin(tmp); -sin(tmp) cos(tmp)]'*dp(k,:)')';
end

% Add noise and drift
dp = dp+sqrt(sp2)*randn(size(dp));
dpsi = dpsi + bias*diff(t(:)) + sqrt(sh2)*randn(length(dpsi),1); 
omega = omega + bias + sqrt(sh2)*randn(length(dpsi),1); 

% Compute position purely from generated odometry
pos_odo = zeros(size(pos_gt,1),2);
pos_odo(1,:) = pos_gt(1,1:2);
for k = 2:size(pos_odo,1)
  tmp = sum(dpsi(1:k-2));
  pos_odo(k,:) = (pos_odo(k-1,:)' + [cos(tmp) sin(tmp); -sin(tmp) cos(tmp)] * dp(k-1,:)')';
end

% Visualise ground truth position and odometry
if visualiseOutput
    figure(1); clf; hold on
    plot(pos_gt(:,1),pos_gt(:,2),'-g')
    plot(pos_odo(:,1),pos_odo(:,2),'-b')
    axis equal
end
%% Construct magnetic field measurements in body frame (as inclination assumed to be known)
m_s = data.m; % Magfield data
m_s = m_s(indCutOff:end,:);
m_b = zeros(size(m_s)); % Magfield data in body frame (so inclination corrected)
for i=1:size(m_s,1) 
    Rincl = roty(-eul(i,2)) * rotx(-eul(i,3));
    m_b(i,:) = Rincl*rotz(-pi/2)*m_s(i,:)'; % Rotate between sensor frame mag and ARkit
end
% Adjust timing because we first do a time update and only then a measurement update
m_b = m_b(2:end,:); 

end
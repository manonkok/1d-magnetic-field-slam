function [MF,MS,PF,GS,MP,PP,undoLoopClosure] = run_filter_from_scratch(m0,P0,Q,dt,dp,omega,loop_start,loop_end,sigma2,gamma_ml,visualiseOutput)

% run_filter_from_scratch - Rerun filter with added loop closure when detected
%
% Syntax:
%   [MF,MS,PF,GS,MP,PP,undoLoopClosure] = run_filter_from_scratch
%       (m0,P0,Q,dt,t,dp,omega,loop_start,loop_end,sigma2,gamma_ml,visualiseOutput)
%
% In:
%   m0              - Initial state estimate
%   P0              - Initial covariance
%   Q               - Process noise covariance
%   dt              - Time between measurements
%   dp              - Delta position measurements
%   omega           - Angular velocity measurements
%   loop_start      - Start indices of detected loops
%   loop_end        - End indices of detected loops
%   loop_end        - End indices of detected loops
%   sigma2          - Covariance of measurement update loop closures
%   gamma_ml        - Threshold to undo erroneous loop closures
%   visualiseOutput - Flag to indicate if outputting messages while running code
%
% Out:
%   MF              - Struct with filtered state estimates
%   MS              - Struct with smoothed state estimates
%   PF              - Struct with filtered state covariances
%   GS              - Struct with cross-covariances 
%   MP              - Struct with predicted state estimates
%   PP              - Struct with predicted state covariances
%   undoLoopClosure - Flag to indicate if outputting messages while running code
%
% Description:
%   Rerun filter with added loop closure when detected. See [1] for details.
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin  

%% Initialisation and pre-allocation
undoLoopClosure = 0; % Flag to detect erroneous loop closures
num_loops = numel(loop_start); % Number of detected loops

% Initialisation EKF state and covariance
% Initial state vector is augmented with landmarks with a prior location 
% of zero with a large uncertainty
m = [m0; zeros(2*num_loops,1)]; 
P = blkdiag(P0,1e4*eye(2*num_loops)); 
  
% Store results as cell arrays
MF = cell(1,size(dp,1));
MS = cell(1,size(dp,1));
PF = cell(1,size(dp,1));
GS = cell(1,size(dp,1));
MP = cell(1,size(dp,1));
PP = cell(1,size(dp,1));
  
% Run EKF until the last loop closure is detected
for k=1:max(loop_end)

    % Previous filter P
    Pf = P;

    % Prediction step
    [m,F,G] = dynamics(m,dp(k,:)',omega(k),dt(k));
    P=F*P*F'+G*Q*G';

    % Make sure the covariance stays PSD
    P=(P+P')/2;    

    % Predicted state mean, covariance, and smoother gain
    mp = m;
    Pp = P;    
    Gs = (Pf*F')/Pp;
    
    % Open loop if any loops to open here
    if any(loop_start==k)
        
        % Make sure not many loops at same point
        if sum(loop_start==k)>1
            error('Opening many loops at same point not supported.');
        end
        
        % Find index
        indLoopClosure = find(loop_start==k,1,'first');

        % Choose elements
        ii = numel(m0)+2*(indLoopClosure-1)+(1:2);
        
        % Measurement model
        H = zeros(2,numel(m));
        H(:,1:2) = eye(2); H(:,ii) = -eye(2);
        R = sigma2*eye(2); 
        
        % Kalman update
        S = H*P*H' + R;
        K = P * H' / S;
        v = ([0; 0] - H * m);
        m = m + K * v;
        P = (eye(size(P))-K*H)*P*(eye(size(P))-K*H)' + K*R*K'; 
    end
    
    % Close loop if any loops to close here
    if any(loop_end==k)
        
        % Make sure not many loops at same point
        if sum(loop_end==k)>1
            error('Closing many loops at same point not supported.');
        end

        % Find index
        indLoopClosure = find(loop_end==k,1,'first');

        % Choose elements
        ii = numel(m0)+2*(indLoopClosure-1)+(1:2);
        
        % Measurement model
        H = zeros(2,numel(m));
        H(:,1:2) = eye(2); H(:,ii) = -eye(2);
        R = sigma2*eye(2);
        
        % Kalman update
        S = H*P*H' + R;
        K = P * H' / S;
        v = ([0; 0] - H * m);
        
        % Check if loop closure is erroneous based on marginal likelihood
        if k == max(loop_end)
            ml = exp(-1/2 * v' / S * v) / (2 * pi * sqrt(det(S)));
            if visualiseOutput
                disp(['ML: ' num2str(ml)])
            end

            if ml < gamma_ml
                if visualiseOutput
                    disp('Undo loop closure')
                end
                undoLoopClosure = 1;
                
                break;
            end
        end
        m = m + K * v;
        P = (eye(size(P))-K*H)*P*(eye(size(P))-K*H)' + K*R*K'; 
    end

    % Store estimate
    MP{k} = mp;
    PP{k} = Pp;
    GS{k} = Gs;
    MF{k} = m;
    PF{k} = P;  
end  
  
  
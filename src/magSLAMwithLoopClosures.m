function [MF,PF,xs,loop_start,loop_end,wp,wm,m_b,t,pos_odo,pos_gt] = ...
    magSLAMwithLoopClosures(filename,driftNoiseParams,makePlots,visualiseOutput,makeVideo)

% magSLAMwithLoopClosures - EKF and RTS smoother for one-dimensional magnetic
% field SLAM with loop closures
%
% Syntax:
%   [MF,PF,xs,loop_start,loop_end,wp,wm,xs,m_b,t,pos_odo,pos_gt] = ...
%     magSLAMwithLoopClosures(filename,driftNoiseParams,makePlots,visualiseOutput,makeVideo)
%
% In:
%   filename            - Filename to load data
%   driftNoiseParams    - Struct indicating what bias and what noise
%                           variances to use for generation of odometry data
%   makePlots           - Flag if code is run to make final plots paper 
%                           (in this case some data sets are not run to the end)
%   visualiseOutput     - Flag to indicate if plotting while running code
%   makeVideo           - Flag to indicate if making a video of the results
%
% Out:
%   MF              - Struct with filtered state estimates
%   PF              - Struct with filtered state covariances
%   xs              - Array with smoothed position and heading states
%   loop_start      - Start indices of detected loops
%   loop_end        - End indices of detected loops
%   wp              - Position weights at final time instance (for plotting)
%   wm              - Magnetic weights at final time instance (for plotting)
%   m_b             - Heading angle derived from ARKit
%   t               - Sampling times
%   pos_gt          - Position derived from ARKit
%   pos_odo         - Position from odometry only
%
% Description:
%   Run EKF and RTS smoother for one-dimensional magnetic field SLAM with 
%   loop closures. See [1] for details.
%
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Preprocessing and settings
% Extract data from file and pre-process
[dp,omega,m_b,t,pos_gt,pos_odo] = prepareData(filename,driftNoiseParams,visualiseOutput);
dt = diff(t);  

% Pre-allocate variables
count_loops = 0; % Loop book keeping
loop_start = []; % Loop book keeping
loop_end = []; % Loop book keeping
MF = cell(1,size(dp,1)); % Filtered mean over time
MS = cell(1,size(dp,1)); % Smoothed mean over time
PF = cell(1,size(dp,1)); % Filtered covariance over time
GS = cell(1,size(dp,1)); % Smoother gain over time
MP = cell(1,size(dp,1)); % Predicted mean over time
PP = cell(1,size(dp,1)); % Predicted mean over time
xs = zeros(3,size(dp,1)); % Smoothed position and heading state estimate

% Settings 
runRTSsmoother = true;
sm = 3; % Magnetometer measurement noise
N_lag = 50; % Do not loop close on the most recent 50 data points
N_lc = 9; % N_lc + 1 is the number of indices checked before deciding on a loop closure
gamma = 0.25;
gamma_ml = 1E-16;
N_dist = 10;
magThreshold = 3;
sigma2 = .01; % Loop-closing measurement noise var
% Process noise
Q = diag([driftNoiseParams.sp2 * ones(1,2), ...
      driftNoiseParams.sh2]);
if contains(filename,'mall') 
    % Slightly higher covariance matrix due to unmodelled errors
    Q = diag([driftNoiseParams.sp2 * ones(1,2), ...
        driftNoiseParams.sh2/10])*10;
end
if contains(filename,'library') 
    % Higher gamma_ml to avoid unwanted loop closures
    gamma_ml = 1E-1;
end

% Initialisation EKF
m0 = [0, 0, 0, 0]'; % Initial state (posx,posy,psi,bias)
P0 = diag([1E-8 1E-8 1e-8 1e-4]); % Initial covariance (posx,posy,psi,bias)
if 3 * sqrt(P0(4,4)) < driftNoiseParams.bias 
    P0(4,4) = driftNoiseParams.bias^2;
end
m = m0;
P = P0;  
  
% Helpers 
indLoopClosure = -inf;
 
if visualiseOutput          
    h = figure(1); clf

    % Adjust figure box size/ratio
    screenSize = get(0, 'ScreenSize');
    aspectRatio = [16 9];
    width = screenSize(3) * 0.5;  % Set the width to half of the screen width
    height = width * (aspectRatio(2) / aspectRatio(1));  % Calculate the height based on the aspect ratio

    % Set the position of the figure [left bottom width height]
    figurePosition = [screenSize(3)/2 - width/2, screenSize(4)/2 - height/2, width, height];
    set(h, 'Position', figurePosition);

    set(h,'Color','w')
    cmap = 1-gray;
    colormap(cmap)
end

%% Run extended Kalman filter
if makeVideo, VideoFrames(length(dp)) = struct('cdata',[],'colormap',[]); end
for k=1:size(dp,1)
    
    % Prediction step
    Pf = P; % Save filtered covariance for later use
    [m,F,G] = dynamics(m,dp(k,:)',omega(k),dt(k)); % Update state
    P=F*P*F'+G*Q*G'; % Update covariance
    P=(P+P')/2; % Make sure the covariance stays positive semi-definite
    mp = m; % Predicted state mean
    Pp = P; % Predicted covariance
    Gs = (Pf*F')/Pp; % Gain for smoother
    
    % Loop closure detection 
    loopClosure = 0;
    wm = zeros(k,2);
    wp = zeros(k,1);
    if k > N_lag + 1
        % Probabilities based on distance
        swp = mean(sqrt(diag(P(1:2,1:2)))); % Extract mean standard deviation 
        d2 = sum((xs(1:2,1:k-1)-m(1:2)).^2,1); % Compute distances squared
        % Compute weights; 8 because: 1/2 from Gaussian distribution, 
        % multiplied by 1/2 because the vector is 2-dimensional,
        % multiplied by 1/2 because we compare two independent sources of information
        wp(1:k-1) = exp(-d2./swp.^2/8); % Weight for distances
        
        % Probabilities based on magnetic field 
        for i = N_lc+1:k-N_lc
            % Compute squared differences between magnetic field measurements 
            % in two directions
            dm2 = sum((m_b(i-N_lc:i,:) - m_b(k-N_lc:k,:)).^2,2);
            dm2_back = sum([m_b(i+N_lc:-1:i,1:2) + m_b(k-N_lc:k,1:2) , m_b(i+N_lc:-1:i,3) - m_b(k-N_lc:k,3)].^2,2);
            % Compute weights; 12 because: 1/2 from Gaussian distribution, 
            % multiplied by 1/3 because the vector is 3-dimensional,
            % multiplied by 1/2 because we compare two independent sources of information
            wm(i,1) = prod(exp(-dm2/12/sm.^2)); 
            wm(i,2) = prod(exp(-dm2_back/12/sm.^2)); 
        end
        
        % Determine index with maximum weight either in fwd or in bwd
        % direction
        [maxProbLoopClosure,indLoopClosure] = max(max(wm(1:k-N_lag,:),[],2).*wp(1:k-N_lag));

        % If probability over threshold gamma, then do loop closure ...
        if maxProbLoopClosure > gamma
            loopClosure = 1;
            % ... unless too close to previous loop closure ... 
            if ismember(indLoopClosure,loop_start) || ...
                any(k-loop_end < N_dist), loopClosure = false; end
            % ... or when magnetic field too close to Earth magnetic field
            if norm(max(m_b(k-N_lc:k,:)) - min(m_b(k-N_lc:k,:))) < magThreshold, loopClosure = false; 
                if visualiseOutput disp('Not enough excitation'); end 
            end
        end
    end

    if visualiseOutput 
        % Plot magnetic field measurements and overlay possible loop
        % closure locations based on position weight
        subplot(221); cla; hold on
            imagesc(t(1:k),ylim,repmat(abs(wp(:))',2,1))
            plot(t(1:k),m_b(1:k,:))   
            xlim([0-eps t(k)]); ylim([-50 20])
            xlabel('Time, t [s]')
            ylabel('Mag [uT]')
            box on, set(gca,'Layer','top')
            clim([0 1])

        % Plot magnetic field and position weights and the overall weights
        subplot(223); cla; hold on
            plot(t(1:k),wp,'-b')    
            plot(t(1:k),wm,'-r')  
            plot(t(1:k),max(wm(1:end,:),[],2).*wp,'-k','LineWidth',2)  
            ylim([0 1])
            if k>1
                xlim([0 t(k)])
            end
            xlabel('Time, t [s]'), ylabel('Weight')
            box on 
    end

    % EKF measurement update
    if loopClosure
        
        % Add to list
        loop_start = [loop_start indLoopClosure];
        loop_end = [loop_end k]; 
        
        % Re-run filter with additional loop closure
        [MF_lc,MS_lc,PF_lc,GS_lc,MP_lc,PP_lc,flag] = ...
          run_filter_from_scratch(m0,P0,Q,dt,dp,omega,loop_start,loop_end,sigma2,gamma_ml,visualiseOutput);        
        
        if ~flag
            % Use versions from rerun filter  
            MF = MF_lc; MS = MS_lc; PF = PF_lc; GS = GS_lc; MP = MP_lc; PP = PP_lc;
          
            % Re-assign
            m = MF{k};
            P = PF{k};
          
            % To visualise current bias estimate    
            if visualiseOutput
                disp(['Current bias estimate: ' num2str(m(4))])
            end

            % Smoother initialization
            MS = MF(1:k);
            ms = MF{k};
            xs = zeros(3,k);
            xs(:,k) = m(1:3);

            % RTS backward recursion
            for j=size(MS,2)-1:-1:1
                ms = MF{j} + GS{j+1} * (ms - MP{j+1});
                MS{j} = ms;
                xs(:,j) = ms(1:3);
            end
          
            % Count the closures
            count_loops = count_loops+1;
            if visualiseOutput
                disp(['Current number of loops: ' num2str(count_loops)])
            end
        else
            % If undoing loop closure, store predicated and filtered estimates, covariances
            loopClosure = 0;
            MP{k} = mp;
            PP{k} = Pp;
            GS{k} = Gs;
            MF{k} = m;
            PF{k} = P; 
            % and undo adding loop closure
            loop_start = loop_start(1:end-1);
            loop_end = loop_end(1:end-1);
            % Keep track of the current best position and orientation state
            xs(:,k) = m(1:3);
        end
    else
        % If no loop closure, store predicated and filtered estimates, covariances
        MP{k} = mp;
        PP{k} = Pp;
        GS{k} = Gs;
        MF{k} = m;
        PF{k} = P; 
        % Keep track of the current best position and orientation state
        xs(:,k) = m(1:3);
    end
    
    % Online plotting
    if visualiseOutput
        subplot(2,2,[2 4]); cla; hold on
            plot(pos_gt(:,1),pos_gt(:,2),'--g')
            plot(xs(1,1:k),xs(2,1:k),'-b')
            if numel(MF{k})>4
                plot(MF{k}(5:2:end),MF{k}(6:2:end),'ob')
            end
            axis square equal
            if loopClosure
                plot(xs(1,indLoopClosure),xs(2,indLoopClosure),'g*')
            end  
            if any(wm(:,1) > 0.7)
                indFwd = find(wm(1:k-N_lag,1) > 0.5);
                plot(xs(1,indFwd),xs(2,indFwd),'b*')
            end
            if any(wm(:,2) > 0.7)
                indBwd = find(wm(1:k-N_lag,2) > 0.5);
                plot(xs(1,indBwd),xs(2,indBwd),'r*')
            end
      
            % Uncertainty radius
            if ~isempty(xs)
                th = linspace(-pi,pi,32);
                uv = diag([sqrt(P(1,1)),sqrt(P(2,2))])*...
                    [cos(th); sin(th)] + xs(1:2,k);
                plot(uv(1,:),uv(2,:),'-b')
        
                % Estimated heading angle (the pi/2 is just for visualization)
                th = pi/2-xs(3,k);
        
                % Heading
                plot(xs(1,k)+[0 cos(th)], xs(2,k)+[0 sin(th)],'-b')
        
            end
            box on
        drawnow
        if makeVideo
            VideoFrames(k) = getframe(gcf);
        end
    end
    % For making the plots in the paper, stop earlier to better show the
    % workings of the algorithm for the square data set
    % if makePlots && strcmp(filename,'../data/short-segments/data-2019-05-07-12-58-29.mat') && (k == 734)
    if makePlots && strcmp(filename,'data/square.mat') && (k == 734)
        VideoFrames = VideoFrames(1:k);
        break;
    end
      
    % For making the plots in the paper, stop earlier to better show the
    % workings of the algorithm for the eight-shaped data set
    % if makePlots && strcmp(filename,'../data/short-segments/data-2019-05-07-13-08-29.mat') && (k == 234)
    if makePlots && strcmp(filename,'data/eight.mat') && (k == 234)    
        VideoFrames = VideoFrames(1:k);
        break;
    end
end
% toc;

if makeVideo
    vidObj = VideoWriter(['1d-magnetic-slam-' filename(6:end-4) '.mp4'],'MPEG-4'); 
    open(vidObj);
    for i=1:length(VideoFrames)
        writeVideo(vidObj,VideoFrames(i));
    end
    close(vidObj);vidObj = VideoWriter(['1d-magnetic-slam-' filename(6:end-4) '.mp4'],'MPEG-4');  
end

end

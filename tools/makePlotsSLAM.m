function [rms_dr,rms_ekf] = makePlotsSLAM(filename,MF,PF,xs,loop_start,loop_end,wp,wm,m_b,t,pos_odo,pos_gt,savePlots)

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

k = length(wm);

%% Plot for square trajectory
if strcmp(filename,'data/square.mat')
    figure(3); clf;
    subplot(221); cla; hold on
    plot(xs(1,:),xs(2,:),'-b')
        if numel(MF{k})>4
            plot(MF{k}(5:2:end),MF{k}(6:2:end),'ob')
        end
        axis equal
        th = linspace(-pi,pi,32);
        uv = chol(PF{k}(1:2,1:2))*[cos(th); sin(th)] + xs(1:2,k);
        plot(uv(1,:),uv(2,:),'-k','LineWidth',2)
        xlabel('x [m]')
        ylabel('y [m]')
        box on 
    subplot(2,2,[3,4]); cla; hold on
      plot(t(1:k),wp,'-b')    
      plot(t(1:k),wm(:,1),'-r')  
      plot(t(1:k),wm(:,2), 'color', [0 0.5 0])  
      plot(t(1:k),max(wm,[],2).*wp,'-k','LineWidth',1.5)  
      ylim([0 1])
      if k>1
        xlim([0 t(k)])
      end
      xlabel('Time [s]'), ylabel('Weights')
      box on 
   subplot(2,2,2); cla; hold on
      plot(t(1:k),wp,'-b')    
      plot(t(1:k),wm(:,1),'-r')  
      plot(t(1:k),wm(:,2), 'color', [0 0.5 0])  
      plot(t(1:k),max(wm,[],2).*wp,'-k','LineWidth',1.5)  
      ylim([0 1])
      if k>1
        xlim([17 21])
      end
      xlabel('Time [s]'), ylabel('Weights')
      box on  

      if savePlots
        matlab2tikz('../paper/figures/square.tex','height','\figureheight', ...
            'width','\figurewidth','extraaxisoptions',[...
            'ylabel style={font=\footnotesize},',...
            'y label style={at={(axis description cs:-0.15,.5)},anchor=south},',...
            'xlabel style={font=\footnotesize},',...
            'yticklabel style={font=\footnotesize},',...
            'xticklabel style={font=\footnotesize},']);
      end
end

%% Plot for eight-shaped trajectory
if strcmp(filename,'data/eight.mat')    
    figure(3); clf;
    subplot(311); cla; hold on
    plot(xs(1,:),xs(2,:),'-b')
        if numel(MF{k})>4
            plot(MF{k}(5:2:end),MF{k}(6:2:end),'ob')
        end
        th = linspace(-pi,pi,32);
        uv = chol(PF{k}(1:2,1:2))*[cos(th); sin(th)] + xs(1:2,k);
        plot(uv(1,:),uv(2,:),'-k','LineWidth',2)
        xlabel('x [m]')
        ylabel('y [m]')
        xlim([-8.1999    4.0076])
        box on
    subplot(3,1,2); cla; hold on
      plot(t(1:k),wp,'-b')    
      plot(t(1:k),wm(:,1),'-r')  
      plot(t(1:k),wm(:,2), 'color', [0 0.5 0])  
      plot(t(1:k),max(wm,[],2).*wp,'-k','LineWidth',1.5)  
      ylim([0 1])
      if k>1
        xlim([0 t(k)])
      end
      xlabel('Time [s]'), ylabel('Weights')
      box on
   subplot(3,1,3); cla; hold on
      plot(t(1:k),m_b(1:k,:))
      xlim([0-eps t(k)]); ylim([-50 20])
      xlabel('Time [s]')
      %ylabel('$y_{\text{m}}^\text{b} [\mu T]$')
      box on
    if savePlots
        matlab2tikz('../paper/figures/eight.tex','height','\figureheight', ...
            'width','\figurewidth','extraaxisoptions',[...
            'ylabel style={font=\footnotesize},',...
            'xlabel style={font=\footnotesize},',...
            'yticklabel style={font=\footnotesize},',...
            'xticklabel style={font=\footnotesize},']);
    end  
end

%% Plot for Kamppi (large scale) data set
if strcmp(filename,'data/mall.mat')
    % For Kamppi data set compute distance travelled
%     sum(sqrt(sum(diff(pos_gt(:,1:2)).^2,2))) % Compute distance per time instance and add all up
    figure(3); clf;
  
    subplot(3,3,1)
        plot(pos_gt(:,1),pos_gt(:,2),'-g')
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        ylim([-30 30])
        yticks([-30,0, 30])
        box on
    subplot(3,3,2)    
        plot(pos_odo(:,1),pos_odo(:,2),'-r')
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        ylim([-30 30])
        yticks([-30,0, 30])
        box on
    subplot(3,3,3), hold on
        plot(xs(1,:),xs(2,:),'-b')
        if numel(MF{k})>4
            plot(MF{k}(5:2:end),MF{k}(6:2:end),'ob')
        end
        xlabel('x [m]')
        ylabel('y [m]')
        ylim([-30 30])
        yticks([-30,0, 30])
        axis equal
        box on
    subplot(3,3,[4 5 6 7 8 9]); hold on 
        for j=1:numel(loop_start)
            x = t(loop_start(j)-[9 0]);
            y = [-50 10];
            fill([x(1) x(2) x(2) x(1) x(1)], ...
                [y(1) y(1) y(2) y(2) y(1)],1,'EdgeColor','w','FaceColor',[.8 .8 .8])
        end
        set(gca,'ColorOrderIndex',1)
        plot(t(1:k),m_b)
        yshift = -110;
        for j=1:numel(loop_start)
            x = t(loop_end(j)-[9 0]);
            y = [-50 10] + yshift;
            fill([x(1) x(2) x(2) x(1) x(1)], ...
                [y(1) y(1) y(2) y(2) y(1)],1,'EdgeColor','w','FaceColor',[.8 .8 .8])
      
            plot([t(loop_start(j)-5) t(loop_end(j)-5)], ...
                [-50 10+yshift],'-','Color',[.8 .8 .8])
        end
        set(gca,'ColorOrderIndex',1)
        plot(t(1:k),m_b+yshift)  
        set(gca,'YTickLabel',[]);
        plot([t(1),t(k)],(mean(m_b(:,3)) + (mean(m_b(:,1)+yshift) - mean(m_b(:,3)))/2)*ones(2,1),'k--')
        xlabel('Time [s]')
        ylabel('$y_{\text{m}}^\text{b}$')
        xlim([0 t(end)])
        box on
    if savePlots
        matlab2tikz('../paper/figures/kamppi.tex','height','\figureheight', ...
            'width','\figurewidth','extraaxisoptions',[...
            'ylabel style={font=\footnotesize},',...
            'y label style={at={(axis description cs:-0.1,.5)},anchor=south},',...
            'x label style={at={(axis description cs:0.5,-.45)},anchor=south},',...
            'xlabel style={font=\footnotesize},',...
            'yticklabel style={font=\footnotesize},',...
            'xticklabel style={font=\footnotesize},']);
    end

    [~, xs_procr] = procrustes(pos_gt(2:end,1:2), xs(1:2,:)');
    [~, pos_odo_procr] = procrustes(pos_gt(:,1:2), pos_odo(:,1:2));
    rms([pos_gt(:,1) - pos_odo_procr(:,1) ; pos_gt(:,2) - pos_odo_procr(:,2)]) % Odometry error
    rms([pos_gt(2:end,1) - xs_procr(:,1) ; pos_gt(2:end,2) - xs_procr(:,2)]) % SLAM error
end

%% Plot for library (small scale) data set
if strcmp(filename,'data/library.mat')    
    figure(3); clf;
  
    subplot(3,3,1)
        plot(pos_gt(:,1),pos_gt(:,2),'-g')
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        box on
    subplot(3,3,2)    
        plot(pos_odo(:,1),pos_odo(:,2),'-r')
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        box on
    subplot(3,3,3), hold on
        plot(xs(1,:),xs(2,:),'-b')
        if numel(MF{k})>4
            plot(MF{k}(5:2:end),MF{k}(6:2:end),'ob')
        end
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        box on
    subplot(3,3,[4 5 6 7 8 9]); hold on 
        for j=1:numel(loop_start)
            x = t(loop_start(j)-[9 0]);
            y = [-50 10];
            fill([x(1) x(2) x(2) x(1) x(1)], ...
                [y(1) y(1) y(2) y(2) y(1)],1,'EdgeColor','w','FaceColor',[.8 .8 .8])
        end
        set(gca,'ColorOrderIndex',1)
        plot(t(1:k),m_b)
        yshift = -100;
        for j=1:numel(loop_start)
            x = t(loop_end(j)-[9 0]);
            y = [-50 10] + yshift;
            fill([x(1) x(2) x(2) x(1) x(1)], ...
                [y(1) y(1) y(2) y(2) y(1)],1,'EdgeColor','w','FaceColor',[.8 .8 .8])
      
            plot([t(loop_start(j)-5) t(loop_end(j)-5)], ...
                [-50 10+yshift],'-','Color',[.8 .8 .8])
        end
        set(gca,'ColorOrderIndex',1)
        plot(t(1:k),m_b+yshift)  
        set(gca,'YTickLabel',[]);
        plot([t(1),t(k)],(mean(m_b(:,3)) + (mean(m_b(:,1)+yshift) - mean(m_b(:,3)))/2)*ones(2,1),'k--')
        xlabel('Time [s]')
        ylabel('Magnetometer readings')
        box on
    if savePlots
        matlab2tikz('../paper/figures/library.tex','height','\figureheight', ...
            'width','\figurewidth','extraaxisoptions',[...
            'ylabel style={font=\footnotesize},',...
            'y label style={at={(axis description cs:-0.12,.5)},anchor=south},',...
            'xlabel style={font=\footnotesize},',...
            'yticklabel style={font=\footnotesize},',...
            'xticklabel style={font=\footnotesize},']);
    end

    % For library data set compute RMSE
    [~, xs_procr] = procrustes(pos_gt(2:end,1:2), xs(1:2,:)');
    [~, pos_odo_procr] = procrustes(pos_gt(:,1:2), pos_odo(:,1:2));
    rms([pos_gt(:,1) - pos_odo_procr(:,1) ; pos_gt(:,2) - pos_odo_procr(:,2)]) % Odometry error
    rms([pos_gt(2:end,1) - xs_procr(:,1) ; pos_gt(2:end,2) - xs_procr(:,2)]) % SLAM error
end

[~, xs_procr] = procrustes(pos_gt(2:k+1,1:2), xs(1:2,1:k)');
[~, pos_odo_procr] = procrustes(pos_gt(:,1:2), pos_odo(:,1:2));
rms_dr = rms([pos_gt(:,1) - pos_odo_procr(:,1) ; pos_gt(:,2) - pos_odo_procr(:,2)])
rms_ekf = rms([pos_gt(2:k+1,1) - xs_procr(:,1) ; pos_gt(2:k+1,2) - xs_procr(:,2)])

end
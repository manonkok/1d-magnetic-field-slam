%% Main file to run one-dimensional magnetic field SLAM algorithm
% References:
%
%   [1] Manon Kok and Arno Solin. Online One-Dimensional Magnetic Field SLAM 
%   with Loop-Closure Detection
%
% Copyright:
%   2024-   Manon Kok and Arno Solin

%% Load data
addpath('tools')
addpath('src')
rng(0,'twister')

clear, clc, close all

% Filenames
filenames = {};
filenames{1} = 'data/library.mat'; 
filenames{2} = 'data/square.mat'; 
filenames{3} = 'data/eight.mat'; 
filenames{4} = 'data/mall.mat'; 

% Settings
makePlots = 1;
makeResults = 1;
runMC = 0;
makePlotsMC = 0;
doComparison = 0;
makeVideo = 0;
visualiseOutput = 1;
savePlots = 0;

%% Run main algorithm
if makeResults
    % Drift and noise settings
    driftNoiseParams.bias = 0.005; %rad/s
    driftNoiseParams.sh2 = 1E-4;
    driftNoiseParams.sp2 = 1E-4; 

    indDataSet = [1;2;3;4];
    rms_drs = [];
    rms_ekfs = [];
    for ind = 1:length(indDataSet)
        [MF,PF,xs,loop_start,loop_end,wp,wm,m_b,t,pos_odo,pos_gt] = ...
            magSLAMwithLoopClosures(filenames{indDataSet(ind)},driftNoiseParams,makePlots,visualiseOutput,makeVideo);
        
        % Visualization
        if makePlots
            [rms_dr,rms_ekf] = makePlotsSLAM(filenames{indDataSet(ind)},MF,PF,xs,loop_start,loop_end,wp,wm,m_b,t,pos_odo,pos_gt,savePlots);
            rms_drs = [rms_drs ; rms_dr];
            rms_ekfs = [rms_ekfs ; rms_ekf];
        end
    end
end

%% Run Monte Carlo simulations
if runMC
    rng(0,'twister')
    indDataSet = 1; % Do Monte Carlo simulations for the library data set
    rms_drs = [];
    rms_ekfs = [];
    
    nMC = 100; % 100 MC simulations per case
    % Consider three cases:
    % 1) Bias different values but stds pos and heading as in paper at 1E-4
    % 2) Different values pos std but bias and heading std as in paper (0.005, 1E-4)
    % 3) Different values heading std but bias and pos std as in paper (0.005, 1E-4)
    biasMC = [0;0.001;0.005;0.01;0.05;0.1]; 
    varPosMC = [1E-6;1E-5;1E-4;1E-3;0.005;1E-2;1E-1];
    varHeadMC = [1E-6;1E-5;1E-4;1E-3;0.005;1E-2;1E-1];
    driftNoiseParamsMCs = cell(3,1);
    driftNoiseParamsMCs{1} = [biasMC,1E-4 * ones(length(biasMC),1), ...
        1E-4 * ones(length(biasMC),1)];
    driftNoiseParamsMCs{2} = [0.005*ones(length(varPosMC),1), ...
        varPosMC, 1E-4 * ones(length(varPosMC),1)];
    driftNoiseParamsMCs{3} = [0.005*ones(length(varHeadMC),1),...
        1E-4 * ones(length(varHeadMC),1),varHeadMC];
    for iFig = 1:length(driftNoiseParamsMCs)
        driftNoiseParamsMC = driftNoiseParamsMCs{iFig};
        rms_ekfs = zeros(size(driftNoiseParamsMC,1),nMC);
        rms_drs = zeros(size(driftNoiseParamsMC,1),nMC);
        for ind = 1:length(indDataSet)
            for iNoise = 1:size(driftNoiseParamsMC,1)
                driftNoiseParams.bias = driftNoiseParamsMC(iNoise,1); %rad/s
                driftNoiseParams.sh2 = driftNoiseParamsMC(iNoise,2);
                driftNoiseParams.sp2 = driftNoiseParamsMC(iNoise,3); 
                for iMC = 1:nMC
                    [MF,PF,xs,loop_start,loop_end,wp,wm,m_b,t,pos_odo,pos_gt] = ...
                        magSLAMwithLoopClosures(filenames{indDataSet(ind)},driftNoiseParams,makePlots,0);
                    [~, xs_procr] = procrustes(pos_gt(2:end,1:2), xs(1:2,:)');
                    [~, pos_odo_procr] = procrustes(pos_gt(:,1:2), pos_odo(:,1:2));
                    rms_dr = rms([pos_gt(:,1) - pos_odo_procr(:,1) ; pos_gt(:,2) - pos_odo_procr(:,2)]);
                    rms_ekf = rms([pos_gt(2:end,1) - xs_procr(:,1) ; pos_gt(2:end,2) - xs_procr(:,2)]);
                    rms_drs(iNoise,iMC) = rms_dr;
                    rms_ekfs(iNoise,iMC) = rms_ekf;
                end
            end
        end
        if iFig == 1
            dataName = 'Bias';
        elseif iFig == 2
            dataName = 'VarPos';
        elseif iFig == 3
            dataName = 'VarHead';
        end
        save(['mcResults' dataName],'rms_ekfs','rms_drs','driftNoiseParamsMC')
    end
end

%% Make plots Monte Carlo simulations with different bias and noise variances
if makePlotsMC
    for iFig = 1:3
        if iFig == 1
            load('mcResultsBias.mat');
            labelString = cell(1,size(driftNoiseParamsMC,1));
            for iNoise = 1:size(driftNoiseParamsMC,1)
                labelString{iNoise} = num2str(driftNoiseParamsMC(iNoise,1));
            end
        elseif iFig == 2
            load('mcResultsVarPos.mat');
            labelString = cell(1,size(driftNoiseParamsMC,1));
            for iNoise = 1:size(driftNoiseParamsMC,1)
                labelString{iNoise} = num2str(driftNoiseParamsMC(iNoise,2));
            end
        elseif iFig == 3
            load('mcResultsVarHead.mat');
            labelString = cell(1,size(driftNoiseParamsMC,1));
            for iNoise = 1:size(driftNoiseParamsMC,1)
                labelString{iNoise} = num2str(driftNoiseParamsMC(iNoise,3));
            end    
        end
        figure; clf;
        boxplot(rms_ekfs','Labels',labelString)
        ax = gca;
        ax.YAxis.Scale ="log";
        ylabel('RMSE [m]')
        if iFig == 1
            xlabel('b_\omega [rad/s]')
            figureName = 'MC_bias';
            xtickLabels = '{0,1E-3,5E-3,1E-2,5E-2,1E-1}';
        elseif iFig == 2
            xlabel('$\sigma_\text{p}^2$ [m$^2$]')
            figureName = 'MC_posVar';
            xtickLabels = '{1E-6,1E-5,1E-4,1E-3,5E-2,1E-2,1E-1}';
        elseif iFig == 3
            xlabel('$\sigma_\omega^2$ [rad$^2$/s$^2$]')
            figureName = 'MC_headingVar';
            xtickLabels = '{1E-6,1E-5,1E-4,1E-3,5E-2,1E-2,1E-1}';
        end
        if savePlots
            figurePath = ['../paper/figures/' figureName '.tex']; 
            matlab2tikz(figurePath,'height','\figureheight', ...
                'width','\figurewidth','extraaxisoptions',[...
                'ylabel style={font=\footnotesize},',...
                'xlabel style={font=\footnotesize},',...
                'yticklabel style={font=\footnotesize},',...
                'xticklabel style={font=\footnotesize},',...
                'xticklabels = ' xtickLabels]);
        end
    end
end
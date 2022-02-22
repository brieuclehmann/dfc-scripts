%% Measurement noise script
% This script simulates data from a fixed dynamic connectivity structure.
% The amplitude of the white noise added to the time series is varied and
% the subsequent effects on the dynamics of functional connectivity are
% analysed. We apply a sliding-window analysis to the simulated fMRI-like 
% time series and use the standard deviation (SD) of the resulting 
% correlation time series as our proxy for the dynamics.

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

rng('default') % for reproducibility

%% Parameter initialisation
% Fix the module structure for each region
toyStruct = [1 1 1 2;1 1 2 2];
nRegions = size(toyStruct,2);
nPairs = (nRegions-1)*nRegions/2;

% Set the number of repetitions/subjects in each group
nReps = 100;

% Generate the brain-state sequence for each subject (same across subjects)
stateSeq = repmat([1;2],1,nReps);

% Set the range of white noise amplitudes
noise = 0:0.2:2.4;

% Initialise the output matrix
sdDFCzcorr = zeros(length(noise),nPairs,nReps);

%% Data generation
for ind = 1:length(noise)
    % Set the amplitude of white noise
    aNoise = noise(ind);
    
    % Generate event sequences for each region
    eventSeq = GenEventSeq(stateSeq,toyStruct);
    
    % Convolve with HRF and add noise to get fMRI-like time series
    TC = GenTC(eventSeq,'aNoise',aNoise);
    
    % Apply sliding window analysis to get pairwise correlation time series
    DFCzcorr = CalcDFCzcorr(TC);
    
    % Calculate standard deviation of correlation time series
    sdDFCzcorr(ind,:,:) = squeeze(std(DFCzcorr,0,2));

end

% Get data for three region-pairs of interest
sdSameMod = squeeze(sdDFCzcorr(:,1,:));
sdDynMod = squeeze(sdDFCzcorr(:,2,:));
sdDiffMod = squeeze(sdDFCzcorr(:,4,:));


%% Statistics
% Prepare data
ySameMod = reshape(sdSameMod',numel(sdSameMod),1);
yDynMod = reshape(sdDynMod',numel(sdSameMod),1);
yDiffMod = reshape(sdDiffMod',numel(sdSameMod),1);

tblSameMod = table(ySameMod,repelem(noise',nReps),'VariableNames',{'Y','Noise'});
tblDynMod = table(yDynMod,repelem(noise',nReps),'VariableNames',{'Y','Noise'});
tblDiffMod = table(yDiffMod,repelem(noise',nReps),'VariableNames',{'Y','Noise'});

% Fit linear models
lmSameMod = fitlm(tblSameMod,'Y ~ Noise^2');
lmDynMod = fitlm(tblDynMod,'Y ~ Noise^2');
lmDiffMod = fitlm(tblDiffMod,'Y ~ Noise^2');

lmCoefs = num2str([lmSameMod.Coefficients.Estimate lmDiffMod.Coefficients.Estimate lmDynMod.Coefficients.Estimate],'%.3f ');
fittedSameMod = unique(lmSameMod.Fitted,'stable');
fittedDiffMod = unique(lmDiffMod.Fitted,'stable');
fittedDynMod = unique(lmDynMod.Fitted,'stable');

%% Plots

gcf = figure;
scatter(repelem(noise,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(noise)]),reshape(sdSameMod',length(noise)*nReps,1),3,'b')
hold on
h1 = plot(noise,fittedSameMod,'-b','LineWidth',1.5);
h2 = scatter(noise,mean(sdSameMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
ax = gca;
set(ax,'xlim',[noise(1)-0.1 noise(end)+0.1])
xlabel('Noise-to-signal ratio (NSR)')
ylabel('SD of correlation \sigma_{corr}')
title('Within-module (R1-R2)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
scatter(repelem(noise,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(noise)]),reshape(sdDiffMod',length(noise)*nReps,1),3,'r')
hold on
h1 = plot(noise,fittedDiffMod,'-r','LineWidth',1.5);
h2 = scatter(noise,mean(sdDiffMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
ax = gca;
set(ax,'xlim',[noise(1)-0.1 noise(end)+0.1])
xlabel('Noise-to-signal ratio (NSR)')
ylabel('SD of correlation \sigma_{corr}')
title('Between-module (R1-R3)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
scatter(repelem(noise,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(noise)]),reshape(sdDynMod',length(noise)*nReps,1),3,'g')
hold on
h1 = plot(noise,fittedDynMod,'-g','LineWidth',1.5);
h2 = scatter(noise,mean(sdDynMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
ax = gca;
set(ax,'xlim',[noise(1)-0.1 noise(end)+0.1])
xlabel('Noise-to-signal ratio (NSR)')
ylabel('SD of correlation \sigma_{corr}')
title('Module change (R1-R4)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow
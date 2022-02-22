%% Neural autocorrelation script
% This script simulates data from a fixed dynamic connectivity structure.
% The autocorrelation of the region-specific event sequences is varied and
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
modStruct = [1 1 1 2;1 1 2 2];
nRegions = size(modStruct,2);
nPairs = (nRegions-1)*nRegions/2;

% Set the number of repetitions/subjects in each group
nReps = 100;

% Generate the brain-state sequence for each subject (same across subjects)
stateSeq = repmat([1;2],1,nReps);

% Set the range of region-specific event sequence autocorrelation
ACrange = -0.8:0.2:0.8;

% Initialise the output matrix
sdDFCzcorr = zeros(length(ACrange),nPairs,nReps);

%% data generation
for ind = 1:length(ACrange)
    % Set the autocorrelation of the region-specific event sequences
    AC = ACrange(ind);
    
    % Generate event sequences for each region
    eventSeq = GenEventSeq(stateSeq,modStruct,'regionAC',AC);
    
    % Convolve with HRF and add noise to get fMRI-like time series
    TC = GenTC(eventSeq,'aNoise',0.2);
    
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

tblSameMod = table(ySameMod,repelem(ACrange',nReps),'VariableNames',{'Y','ACrange'});
tblDynMod = table(yDynMod,repelem(ACrange',nReps),'VariableNames',{'Y','ACrange'});
tblDiffMod = table(yDiffMod,repelem(ACrange',nReps),'VariableNames',{'Y','ACrange'});

% Fit linear models
lmSameMod = fitlm(tblSameMod,'Y ~ ACrange^2');
lmDynMod = fitlm(tblDynMod,'Y ~ ACrange^2');
lmDiffMod = fitlm(tblDiffMod,'Y ~ ACrange^2');

lmCoefs = num2str([lmSameMod.Coefficients.Estimate lmDiffMod.Coefficients.Estimate lmDynMod.Coefficients.Estimate],'%.3f ');
fittedSameMod = unique(lmSameMod.Fitted,'stable');
fittedDiffMod = unique(lmDiffMod.Fitted,'stable');
fittedDynMod = unique(lmDynMod.Fitted,'stable');
lmSameModCI = coefCI(lmSameMod);
lmDiffModCI = coefCI(lmDiffMod);
lmDynModCI = coefCI(lmDynMod);

%% Plots

gcf = figure;
scatter(repelem(ACrange,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(ACrange)]),reshape(sdSameMod',length(ACrange)*nReps,1),3,'b')
hold on
h1 = plot(ACrange,fittedSameMod,'-b','LineWidth',1.5);
h2 = scatter(ACrange,mean(sdSameMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
xlabel('Event autocorrelation \rho_{reg}')
ylabel('SD of correlation \sigma_{corr}')
title('Within-module (R1-R2)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
scatter(repelem(ACrange,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(ACrange)]),reshape(sdDiffMod',length(ACrange)*nReps,1),3,'r')
hold on
h1 = plot(ACrange,fittedDiffMod,'-r','LineWidth',1.5);
h2 = scatter(ACrange,mean(sdDiffMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
xlabel('Event autocorrelation \rho_{reg}')
ylabel('SD of correlation \sigma_{corr}')
title('Between-module (R1-R3)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
scatter(repelem(ACrange,nReps)+repmat(linspace(-0.02,0.02,nReps),[1,length(ACrange)]),reshape(sdDynMod',length(ACrange)*nReps,1),3,'g')
hold on
h1 = plot(ACrange,fittedDynMod,'-g','LineWidth',1.5);
h2 = scatter(ACrange,mean(sdDynMod,2),'s','LineWidth',1.5,'MarkerEdgeColor','k');
legend([h1 h2],{'Fitted values','Mean SD'},'Location','best')
xlabel('Event autocorrelation \rho_{reg}')
ylabel('SD of correlation \sigma_{corr}')
title('Module change (R1-R4)')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow
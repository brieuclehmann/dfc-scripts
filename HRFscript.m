%% HRF simulation script
% This script simulates data from a fixed dynamic connectivity structure.
% The parameters of the haemodynamic response function (HRF) are varied and
% the subsequent effects on the dynamics of functional connectivity are
% analysed. Here, we vary the delay of the peak response and the dispersion
% of the peak response. We apply a sliding-window analysis to the simulated
% fMRI-like time series and use the standard deviation (SD) of the
% resulting correlation time series as our proxy for the dynamics.

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

% Set the range of HRF parameters
delay_response = 5:0.2:9;
disp_response = 0.6:0.2:2.4;

% Initialise the output matrices
sdDFCzcorr = zeros(length(delay_response),length(disp_response),nPairs,nReps);

%% Data simulation
for ind1 = 1:length(delay_response)
    for ind2 = 1:length(disp_response)
        % Set the HRF parameters
        P1 = delay_response(ind1);
        P3 = disp_response(ind2);
        
        P(1) = P1;     % delay of response (relative to onset)
        P(2) = 15;     % delay of undershoot (relative to onset)
        P(3) = P3;     % dispersion of response
        P(4) = 1;      % dispersion of undershoot
        P(5) = 3;      % ratio of response to undershoot
        P(6) = 0;      % onset (seconds)
        P(7) = 32;     % length of kernel (seconds)
        
        % Generate event sequences for each region
        eventSeq = GenEventSeq(stateSeq,toyStruct);
        
        % Convolve with HRF and add noise to get fMRI-like time series
        TC = GenTC(eventSeq,'P',P);

        % Apply sliding window analysis to get pairwise correlation time series
        DFCzcorr = CalcDFCzcorr(TC);
        
        % Calculate standard deviation of correlation time series
        sdDFCzcorr(ind1,ind2,:,:) = squeeze(std(DFCzcorr,0,2));
    end
end

% Get data for three region-pairs of interest
sdSameMod = squeeze(sdDFCzcorr(:,:,1,:));
sdDynMod = squeeze(sdDFCzcorr(:,:,2,:));
sdDiffMod = squeeze(sdDFCzcorr(:,:,4,:));

%% Statistics
% Prepare data
ySameMod = reshape(permute(sdSameMod,[3 1 2]),numel(sdSameMod),1);
yDynMod = reshape(permute(sdDynMod,[3 1 2]),numel(sdDynMod),1);
yDiffMod = reshape(permute(sdDiffMod,[3 1 2]),numel(sdDiffMod),1);

tblSameMod = table(ySameMod,repelem(delay_response',nReps*length(disp_response)),repelem(repmat(disp_response',[length(delay_response) 1]),nReps),'VariableNames',{'Y','Delay_response','Disp_response'});
tblDynMod = table(yDynMod,repelem(delay_response',nReps*length(disp_response)),repelem(repmat(disp_response',[length(delay_response) 1]),nReps),'VariableNames',{'Y','Delay_response','Disp_response'});
tblDiffMod = table(yDiffMod,repelem(delay_response',nReps*length(disp_response)),repelem(repmat(disp_response',[length(delay_response) 1]),nReps),'VariableNames',{'Y','Delay_response','Disp_response'});

lmSameMod = fitlm(tblSameMod,'Y ~ Delay_response * Disp_response + Delay_response^2 + Disp_response^2');
lmDynMod = fitlm(tblDynMod,'Y ~ Delay_response * Disp_response + Delay_response^2 + Disp_response^2');
lmDiffMod = fitlm(tblDiffMod,'Y ~ Delay_response * Disp_response + Delay_response^2 + Disp_response^2 ');

lmCoefs = num2str([lmSameMod.Coefficients.Estimate lmDiffMod.Coefficients.Estimate lmDynMod.Coefficients.Estimate],'%.3f ');
lmSameModCI = coefCI(lmSameMod);
lmDiffModCI = coefCI(lmDiffMod);
lmDynModCI = coefCI(lmDynMod);

%% Plot results
% Get data for bar plots
sdSameModYoung = squeeze(sdDFCzcorr(delay_response==6,disp_response==1,1,:))';
sdDynModYoung = squeeze(sdDFCzcorr(delay_response==6,disp_response==1,2,:))';
sdDiffModYoung = squeeze(sdDFCzcorr(delay_response==6,disp_response==1,4,:))';

sdSameModOld = squeeze(sdDFCzcorr(delay_response==8,disp_response==2,1,:))';
sdDynModOld = squeeze(sdDFCzcorr(delay_response==8,disp_response==2,2,:))';
sdDiffModOld = squeeze(sdDFCzcorr(delay_response==8,disp_response==2,4,:))';

gcf = figure;
imagesc(disp_response,delay_response,mean(sdSameMod,3)); c = colorbar;
ylabel(c,'Mean SD of correlation'); title('Within module')
ylabel('Delay of response \tau_{HRF}'); xlabel('Dispersion of response \sigma_{HRF}')
set(gca,'YDir','normal'); hold on; plot([1 1],ylim,'--k'); 
plot([2 2],ylim,'--k'); plot(xlim,[6 6],'--k'); plot(xlim,[8 8],'--k')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
imagesc(disp_response,delay_response,mean(sdDiffMod,3)); c = colorbar;
ylabel(c,'Mean SD of correlation'); title('Between module') 
ylabel('Delay of response \tau_{HRF}'); xlabel('Dispersion of response \sigma_{HRF}')
set(gca,'YDir','normal'); hold on; plot([1 1],ylim,'--k')
plot([2 2],ylim,'--k'); plot(xlim,[6 6],'--k'); plot(xlim,[8 8],'--k')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
imagesc(disp_response,delay_response,mean(sdDynMod,3)); c = colorbar;
ylabel(c,'Mean SD of correlation'); title('Module change')
ylabel('Delay of response \tau_{HRF}'); xlabel('Dispersion of response \sigma_{HRF}')
set(gca,'YDir','normal'); hold on; plot([1 1],ylim,'--k')
plot([2 2],ylim,'--k'); plot(xlim,[6 6],'--k'); plot(xlim,[8 8],'--k')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
pos_same = [1 1.2];
pos_diff = [1.5 1.7];
pos_dyn = [2.0 2.2];
box_same = boxplot([sdSameModYoung; sdSameModOld; sdDiffModYoung; sdDiffModOld; sdDynModYoung; sdDynModOld]','colors',[0.494 0.184 0.556; 0.929 0.694 0.125],'positions',[pos_same pos_diff pos_dyn],'width',0.15);
set(box_same,'linewidth',1.5)
hold on
set(gca,'xlim',[0.8 2.4])
set(gca,'XTick',[1.1 1.6 2.1])
set(gca,'XTickLabel',{'R1-R2','R1-R3','R1-R4'})
obj = findall(gca,'Tag','Box');
obj = obj([2;1],:);
hLegend = legend(obj, {'G1','G2'},'Location','best');
plot([1.35 1.35],ylim,'--k')
plot([1.85 1.85],ylim,'--k')
ylabel('SD of correlation')
title('G1: \sigma_{HRF} = 1, \tau_{HRF} = 6      G2: \sigma_{HRF} = 2, \tau_{HRF} = 8');
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.5 2.8];
drawnow
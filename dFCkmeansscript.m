%% k-means analysis script 
% This script simulates data from two groups of individuals: the first
% group can visit 9 FC states while the second group can visit only 6 of
% these 9 FC states. Each individual experiences 3 FC state transitions at 
% regular intervals. We apply a sliding-window analysis to the simulated 
% fMRI-like time series and then use a k-means analysis to try to recover
% the FC states for each individual at each window.

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

rng('default') % for reproducibility

%% Parameter initialisation
% modStruct (9 states x 32 regions) contains the module assignments 
% for each region in each state
load modStruct32
modStruct = modStruct32;
% Alternatively, use GenModStruct to generate different module structures

% Set the number of FC states attainable by individuals in each group
nStatesG1 = size(modStruct,1);
nStatesG2 = 6;

% Set the number of individuals in each group
nSub = 5;
G1Ind = 1:nSub;
G2Ind = nSub+1:nSub+nSub;

% Set the FC state transition matrix for G1, specifiying the probability of
% transitioning from one state to another. Individuals in G2 have an equal
% probability of transitioning from one state to any other state.
transitionMatG1(1:6,:) = repmat([1/10 1/10 1/10 1/10 1/10 1/10 1/6 1/6 1/6],6,1);
transitionMatG1(7:9,:) = repmat([1/12 1/12 1/12 1/12 1/12 1/12 1/4 1/4 1/4],3,1);
% Transitions from one state to itself are not permitted.
transitionMatG1(eye(nStatesG1) ~= 0) = 0;

%% Data generation
% Generate a FC state sequence for each subject in both groups
stateSeqG1 = GenStateSeq(nStatesG1,'nSub',nSub,'transitionMat',transitionMatG1);
stateSeqG2 = GenStateSeq(nStatesG2,'nSub',nSub);
stateSeqAll = cat(2,stateSeqG1,stateSeqG2);

% Generate the event sequences for each subject
eventSeqG1 = GenEventSeq(stateSeqG1, modStruct);
eventSeqG2 = GenEventSeq(stateSeqG2, modStruct);

% Convolve with HRF and add noise to get fMRI-like time series
TCG1 = GenTC(eventSeqG1);
TCG2 = GenTC(eventSeqG2);
TCAll = cat(3,TCG1,TCG2);

% Apply sliding window analysis to get correlation matrices across time
DFCzcorrG1 = CalcDFCzcorr(TCG1);
DFCzcorrG2 = CalcDFCzcorr(TCG2);
DFCzcorrAll = cat(3,DFCzcorrG1,DFCzcorrG2);

% Calculate 'true' (average) correlation matrices for each FC state
avDFCzcorrG1 = trueDFCzcorr(TCG1,stateSeqG1);
avDFCzcorrG2 = trueDFCzcorr(TCG2,stateSeqG2);
avDFCzcorrAll = trueDFCzcorr(TCAll,stateSeqAll);


%% k-means analysis 
% Run k-means clustering for the two groups, and the total combined group
% Combined group
[IDXAll, centroidsAll, sumAll] = kMeansFMRI(DFCzcorrAll);
% Group 1 (9 FC states)
[IDXG1, centroidsG1, sumG1] = kMeansFMRI(DFCzcorrG1);
% Group 2 (6 of these 9 FC states)
[IDXG2, centroidsG2, sumG2] = kMeansFMRI(DFCzcorrG2);

% Compute error measures for the two groups, and the total combined group
% (warning: errorkMeans assumes number of time points nT = 360)
[err_DFCzcorr, err_nDiffStates] = errorkMeans(stateSeqAll, avDFCzcorrAll, IDXAll, centroidsAll);       
[err_DFCzcorrG2, err_nDiffStatesG2] = errorkMeans(stateSeqG2, avDFCzcorrG2, IDXG2, centroidsG2);       
[err_DFCzcorrG1, err_nDiffStatesG1] = errorkMeans(stateSeqG1, avDFCzcorrG1, IDXG1, centroidsG1);

%% Statistics
for k1 = 1:12
    % Test if error in number of FC state transitions is different from
    % zero for G1 in combined analysis
    [~, r_nStatesG1Comb(k1)] = signrank(err_nDiffStates(G1Ind,k1));
    % Test if error in number of FC state transitions is different from
    % zero for G2 in combined analysis
    [~, r_nStatesG2Comb(k1)] = signrank(err_nDiffStates(G2Ind,k1));
    
    % Test if number of FC state transitions is different between G1 and
    % G2 in combined analysis
    [~, r_diffComb(k1)] = ranksum(err_nDiffStates(G2Ind,k1),err_nDiffStates(G1Ind,k1));
    
    % Test if centroid errors are different between G1 and G2 in combined
    % analyses
    [~, r_errorComb(k1)] = ranksum(err_DFCzcorr(G2Ind,k1),err_DFCzcorr(G1Ind,k1));
    % Test if centroid errors are different between G1 and G2 in separate
    % analyses
    [~, r_errorSep(k1)] = ranksum(err_DFCzcorrG2(:,k1),err_DFCzcorrG1(:,k1));
    
    % Test if error in number of FC state transitions is different from
    % zero for G1 in separate analysis
    [~, r_nStatesG1Sep(k1)] = signrank(err_nDiffStatesG1(:,k1));
    
    % Test if error in number of FC state transitions is different from
    % zero for G2 in separate analysis
    [~, r_nStatesG2Sep(k1)] = signrank(err_nDiffStatesG2(:,k1));
    for k2 = 1:12
        % Test if centroid errors are different between G1 and G2 in
        % separate analyses
        [~, r_diffSep(k1,k2)] = ranksum(err_nDiffStatesG2(:,k1),err_nDiffStatesG1(:,k2));
        
        % For plotting later
        diff_nDiffStates(k1,k2) = mean(err_nDiffStatesG2(:,k1))-mean(err_nDiffStatesG1(:,k2));
    end
end

% Prepare data for plotting
r_nStatesG1Comb = double(r_nStatesG1Comb);
r_nStatesG2Comb = double(r_nStatesG2Comb);
r_diffComb = double(r_diffComb);
r_errorComb = double(r_errorComb);
r_errorSep = double(r_errorSep);

r_nStatesG1Comb(r_nStatesG1Comb==0)=NaN;
r_nStatesG2Comb(r_nStatesG2Comb==0)=NaN;
r_diffComb(r_diffComb==0)=NaN;
r_errorComb(r_errorComb==0)=NaN;
r_errorSep(r_errorSep==0)=NaN;

r_nStatesG1Sep = double(r_nStatesG1Sep);
r_nStatesG2Sep = double(r_nStatesG2Sep);
r_diffSep = double(r_diffSep);

r_nStatesG1Sep(r_nStatesG1Sep==0)=NaN;
r_nStatesG2Sep(r_nStatesG2Sep==0)=NaN;
r_diffSep(r_diffSep==0)=NaN;
r_diffSep = bsxfun(@times,(1:12)',r_diffSep);

%% Plots

gcf = figure;
hold on; 
plot(mean(err_nDiffStates(G1Ind,:),1),'Color',[0.929 0.694 0.125],'LineWidth',2);
plot(mean(err_nDiffStates(G2Ind,:),1),'Color',[0.494 0.184 0.556],'LineWidth',2);
scatter(1:12,r_nStatesG2Comb.*mean(err_nDiffStates(G2Ind,:),1),20,[0.494 0.184 0.556],'*','LineWidth',1);
scatter(1:12,r_nStatesG1Comb.*mean(err_nDiffStates(G1Ind,:),1),20,[0.929 0.694 0.125],'*','LineWidth',1);
xlabel('Number of clusters (k)'); ylabel({'Mean error in number', 'of state changes'});
legend('G1 (combined)','G2 (combined)','Location','NorthWest')
set(gca,'xlim',[1 12],'ylim',[-5 10])
line(xlim,[0 0],'Color','k')
title('Combined analysis')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
hold on;
plot(mean(err_nDiffStatesG1,1),'Color',[0.929 0.694 0.125],'LineWidth',2,'LineStyle',':');
plot(mean(err_nDiffStatesG2,1),'Color',[0.494 0.184 0.556],'LineWidth',2,'LineStyle',':');
scatter(1:12,r_nStatesG1Sep.*mean(err_nDiffStatesG1,1),20,[0.929 0.694 0.125],'*','LineWidth',1);
scatter(1:12,r_nStatesG2Sep.*mean(err_nDiffStatesG2,1),20,[0.494 0.184 0.556],'*','LineWidth',1);
xlabel('Number of clusters (k)'); ylabel({'Mean error in number', 'of state changes'});
legend('G1 (separate)','G2 (separate)','Location','NorthWest')
set(gca,'xlim',[1 12],'ylim',[-5 10])
line(xlim,[0 0],'Color','k')
title('Separate analysis')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
plot(mean(err_nDiffStates(G2Ind,:),1)-mean(err_nDiffStates(G1Ind,:),1),'LineWidth',2,'Color','k');
hold on
line(xlim,[0 0],'Color','k')
scatter(1:12,r_diffComb.*(mean(err_nDiffStates(G2Ind,:),1)-mean(err_nDiffStates(G1Ind,:),1)),20,'k','*','LineWidth',1);
set(gca,'xlim',[1 12])
xlabel('Number of clusters (k)');
ylabel({'Group difference in mean number','of state changes: G2 - G1'})
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
axmax = max(max(diff_nDiffStates));
axmin = min(min(diff_nDiffStates));
imagesc([1 12], [1 12], diff_nDiffStates); caxis([-axmax axmax]);
hold on
scatter(repelem(1:12,12),r_diffSep(:),20,'k','*','LineWidth',1);
colorbar; colormap(jet(1000)); set(colorbar,'Limits',[axmin axmax]);
set(gca,'YDir','normal'); xlabel('Number of clusters (k): G1');
ylabel('Number of clusters (k): G2');
title({'Group difference in mean number','of state changes: G2 - G1'})
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
hold on;
plot(mean(err_DFCzcorr(G1Ind,:),1),'Color',[0.929 0.694 0.125],'LineWidth',2);
plot(mean(err_DFCzcorr(G2Ind,:),1),'Color',[0.494 0.184 0.556],'LineWidth',2);
scatter(1:12,r_errorComb.*(mean(err_DFCzcorr(G2Ind,:),1)+mean(err_DFCzcorr(G1Ind,:),1))/2,20,'k','*','LineWidth',1);
xlabel('Number of clusters (k)'); ylabel('Mean error of centroids');
legend('G1 (combined)','G2 (combined)','Location','NorthEast')
set(gca,'xlim',[1 12])
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
hold on;
plot(mean(err_DFCzcorrG1,1),'Color',[0.929 0.694 0.125],'LineWidth',2,'LineStyle',':');
plot(mean(err_DFCzcorrG2,1),'Color',[0.494 0.184 0.556],'LineWidth',2,'LineStyle',':');
scatter(1:12,r_errorSep.*(mean(err_DFCzcorrG2,1)+mean(err_DFCzcorrG1,1))/2,20,'k','*','LineWidth',1);
xlabel('Number of clusters (k)'); ylabel('Mean error of centroids');
legend('G1 (separate)','G2 (separate)','Location','NorthEast')
set(gca,'xlim',[1 12])
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow

gcf = figure;
plot(sumAll/2,'k','LineWidth',2)
hold on
plot(sumG1,'Color',[0.929 0.694 0.125],'LineWidth',2,'LineStyle',':')
plot(sumG2,'Color',[0.494 0.184 0.556],'LineWidth',2,'LineStyle',':')
set(gca,'xlim',[1 12])
xlabel('Number of clusters (k)'); ylabel('Residual sum of squares');
legend('Combined','G1 (separate)','G2 (separate)','Location','NorthEast')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 3.2 2.5];
drawnow
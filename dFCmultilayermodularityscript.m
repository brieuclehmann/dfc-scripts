%% Multilayer modularity analysis script
% This script simulates data from two groups of individuals: those in the 
% first group have one FC state transition while those in the second group
% have 3 FC state transitions. Every individual can visit the same set of 
% 9 FC states. We apply a sliding-window analysis to the simulated 
% fMRI-like time series and then use a multilayer modularity algorithm to 
% try to recover the module structure for each individual at each window.

% Dependency: GenLouvain package - freely available on 
% http://netwiki.amath.unc.edu/GenLouvain/GenLouvain

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

rng('default') % for reproducibility

%% Parameter initialisation
% modStruct (9 states x 32 regions) contains the module assignments 
% for each region in each state
load modStruct32
modStruct = modStruct32;
nRegions = size(modStruct,2);
% Alternatively, use GenModStruct to generate different module structures

% Get the number of FC states attainable by individuals in each group
nStates = size(modStruct,1);

% Set the number of different states for each group - this is the number of
% FC state transitions + 1
nDiffStatesG1 = 2;
nDiffStatesG2 = 4;

% Set the number of individuals in each group
nSub = 50;
G1Ind = 1:nSub;
G2Ind = nSub+1:nSub+nSub;

%% Data generation
% Generate a FC state sequence for each subject in both groups
stateSeqG1 = GenStateSeq(nStates,'nDiffStates',nDiffStatesG1,'nSub',nSub); 
stateSeqG2 = GenStateSeq(nStates,'nDiffStates',nDiffStatesG2,'nSub',nSub);
stateSeqAll = cat(2,repelem(stateSeqG1,2,1),stateSeqG2);

% Generate the event sequences for each subject
eventSeq = GenEventSeq(stateSeqAll, modStruct);

% Convolve with HRF and add noise to get fMRI-like time series
TC = GenTC(eventSeq);

% Apply (non-overlapping) sliding window analysis to get correlation
% matrices across time
winP = [30 30];
DFCzcorr = CalcDFCzcorr(TC,winP);

% Get true module structure for each window
nWin = size(DFCzcorr,2);
nWinState = nWin/size(stateSeqAll,1);
winStateSeq = repelem(stateSeqAll,nWinState,1);
winModStruct = reshape(modStruct(winStateSeq,:)',nRegions,nWin,nSub+nSub);

%% Multilayer modularity analysis
% Set hyperparameter range
Gamma = 0.75:0.25:2.5;
Omega = 0.25:0.25:4;

% Preallcoate output
err_structG1 = zeros(length(Gamma),length(Omega));
err_structG2 = zeros(length(Gamma),length(Omega));
err_flexG1 = zeros(length(Gamma),length(Omega));
err_flexG2 = zeros(length(Gamma),length(Omega));
mmFlexG1 = zeros(length(Gamma),length(Omega));
mmFlexG2 = zeros(length(Gamma),length(Omega));

for g = 1:length(Gamma)
    for w = 1:length(Omega)
        % Run multilayer modularity algorithm
        multiModStruct = multiModFMRI(DFCzcorr,Gamma(g),Omega(w));
        
        % Compute errors for each group
        [err_struct, err_flex, mmFlex, trueFlex] = errorMultiMod(multiModStruct, winModStruct);
        % Get mean error in connectivity structure
        err_structG1(g,w) = mean(err_struct(G1Ind));
        err_structG2(g,w) = mean(err_struct(G2Ind));
        % Get mean error in flexibility
        err_flexG1(g,w) = mean(err_flex(G1Ind));
        err_flexG2(g,w) = mean(err_flex(G2Ind));
        % Get mean recovered flexibility
        mmFlexG1(g,w) = mean(mmFlex(G1Ind));
        mmFlexG2(g,w) = mean(mmFlex(G2Ind));
    end
end
% Get 'true' (average) flexibility for both groups
trueFlexG1 = mean(trueFlex(G1Ind));
trueFlexG2 = mean(trueFlex(G2Ind));
trueFlexDiff = trueFlexG1 - trueFlexG2;

%% Plots
x = [Omega(1) Omega(end)];
y = [Gamma(1) Gamma(end)];

[G1min(1), G1min(2)] = find(err_structG1 == min(err_structG1(:)));
[G2min(1), G2min(2)] = find(err_structG2 == min(err_structG2(:)));
G1min(1) = Gamma(G1min(1));
G1min(2) = Omega(G1min(2));
G2min(1) = Gamma(G2min(1));
G2min(2) = Omega(G2min(2));

axmin = min(min([err_structG1, err_structG2]));
axmax = max(max([err_structG1, err_structG2]));
axstruct = [axmin axmax];

gcf = figure;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.25 .1 .65 .8]);
imagesc(x,y,err_structG1); 
set(gca,'YDir','normal'); xlabel('omega \omega'); ylabel('gamma \gamma');
c = colorbar; drawnow; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; drawnow; ax.Position = axpos;
ylabel(c,{'Mean percentage error in','connectivity structure'})
caxis(axstruct); colormap(jet)
ax.ActivePositionProperty = 'OuterPosition';
hold on
crossG1 = plot(G1min(2),G1min(1),'rx','MarkerSize',10,'LineWidth',3);
legend(crossG1,'Minimum error')
title('G1: error in connectivity structure')
gcf.PaperUnits = 'inches'; gcf.PaperPosition = [0 0 4 2.8];
drawnow

gcf = figure;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.25 .1 .65 .8]);
imagesc(x,y,err_structG2); 
c = colorbar; drawnow; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; drawnow; ax.Position = axpos;
ax.ActivePositionProperty = 'OuterPosition';
ylabel(c,{'Mean percentage error in','connectivity structure'})
caxis(axstruct); colormap(jet)
set(gca,'YDir','normal'); xlabel('omega \omega'); ylabel('gamma \gamma');
hold on
crossG2 = plot(G2min(2),G2min(1),'rx','MarkerSize',10,'LineWidth',3);
legend(crossG2,'Minimum error')
text(Omega(1) - 1.7, 0.5*(Gamma(1)+Gamma(end)),{'G2:','3 state','transitions'},'HorizontalAlignment','center','FontWeight','bold')
title('G2: error in connectivity structure')
gcf.PaperUnits = 'inches';
gcf.PaperPosition = [0 0 4 2.8];
drawnow

axmin = min(min([err_flexG1, err_flexG2]));
axmax = max(max([err_flexG1, err_flexG2]));
axstruct = [-axmax axmax];

gcf = figure;
imagesc(x,y,err_flexG1);
c = colorbar; set(c,'Limits',[axmin axmax]);
drawnow; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; drawnow; ax.Position = axpos;
ax.ActivePositionProperty = 'OuterPosition';
ylabel(c,'Mean error in mean flexibility')
set(gca,'YDir','normal'); xlabel('omega \omega'); ylabel('gamma \gamma');
caxis(axstruct); colormap(jet); title('G1: error in flexibility')
gcf.PaperUnits = 'inches'; gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

gcf = figure;
imagesc(x,y,err_flexG2);
c = colorbar; ylabel(c,'Mean error in mean flexibility')
drawnow; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; drawnow; ax.Position = axpos;
ax.ActivePositionProperty = 'OuterPosition';
set(c,'Limits',[axmin axmax]);
set(gca,'YDir','normal'); xlabel('omega \omega'); ylabel('gamma \gamma');
caxis(axstruct); colormap(jet); title('G2: error in flexibility')
gcf.PaperUnits = 'inches'; gcf.PaperPosition = [0 0 3.5 2.8];
drawnow

axmin = min(min(mmFlexG2-mmFlexG1));
axmax = max(max(mmFlexG2-mmFlexG1));
axstruct = [-axmax axmax];

gcf = figure;
imagesc(x,y,mmFlexG2-mmFlexG1); caxis(axstruct); colormap(jet)
c = colorbar; ylabel(c,'Group difference in mean flexibility (G2 - G1)')
drawnow; ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(3) = 0.5*cpos(3); c.Position = cpos; drawnow; ax.Position = axpos;
ax.ActivePositionProperty = 'OuterPosition';
set(gca,'YDir','normal'); xlabel('omega \omega'); ylabel('gamma \gamma');
set(c,'Limits',[axmin axmax]);
drawnow
gcf.PaperUnits = 'inches'; gcf.PaperPosition = [0 0 3.5 3.2];
drawnow
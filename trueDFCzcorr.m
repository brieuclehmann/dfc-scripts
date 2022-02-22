function avDFCzcorr = trueDFCzcorr(TC,stateSeq)
% Calculate 'true' (average) Fisher-transformed correlation matrices for
% each FC state. This is done by finding each time period when an
% individual is in a given FC state, computing the correlation between all
% region pairs in these time periods, and averaging over all such periods.

% avDFCzcorr is a nStates x nPairs matrix with (ij)th entry 
% corresponding to the average Fisher-transformed correlation between 
% region pair j in FC state i

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
nRegions = size(TC,1);
nPairs = nRegions*(nRegions-1)/2;
nT = size(TC,2);
nSub = size(TC,3);
nDiffStates = size(stateSeq,1);
nStates = max(max(stateSeq));
avDFCzcorr = zeros(nStates, nPairs);

% Get time between FC state transitions
nTState = nT/nDiffStates;
winP = [nTState nTState];

%% Calculate average correlation matrices for each FC state
% Calculate correlation matrices for each state traversed
stateDFCzcorr=CalcDFCzcorr(TC,winP);

dat = reshape(stateDFCzcorr,nPairs,nDiffStates*nSub);

% Compute average correlation matrices across each FC state
for s1 = 1:nStates
    avDFCzcorr(s1,:) = mean(dat(:,stateSeq==s1),2);
end
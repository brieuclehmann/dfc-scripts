function [err_DFCzcorr, err_nDiffStates] = errorkMeans(stateSeq, avDFCzcorr, IDX, centroids)
% Compute error measures for k-means analysis for each k. err_DFCzcorr is 
% the centroid error for each subject while err_nDiffStates is the error in
% number of state transitions for each subject. 

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
nT = 360; % default - careful!
nClust = size(IDX,1);
nSub = size(IDX,2);
nWin = size(IDX,3);

err_DFCzcorr = zeros(nSub,nClust);
err_nDiffStates = zeros(nSub,nClust);

%% Computation of error measures
for n = 1:nSub
    % Get state sequence for this subject
    nDiffStates = sum(diff([stateSeq(:,n)' eps])~=0);
    tState = nT/size(stateSeq,1);
    overlap = (nT - nWin + 1)/2;
    thisSeq = repelem(stateSeq(:,n),tState);
    thisSeq = thisSeq(overlap:nT-overlap);
    for k = 1:nClust
        % Get k-means state sequence for this k
        thisIDX = squeeze(IDX(k,n,:));
        % Compute error in number of state changes
        err_nDiffStates(n,k) = sum(diff([thisIDX' eps])~=0)-nDiffStates;
        
        % Compute error for centroid relative to 'true' correlation matrix
        % for each window
        errorMat = avDFCzcorr(thisSeq,:) - squeeze(centroids(k,thisIDX,:));
        err_DFCzcorr(n, k) = sum(sqrt(sum(abs(errorMat).^2,2)));
    end
end
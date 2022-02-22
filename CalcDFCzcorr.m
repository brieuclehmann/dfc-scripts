function vecDFCzcorr = CalcDFCzcorr(TC, winP)
% Calculate Fisher-transformed sliding window correlations of the BOLD time
% series TC for each pair of regions in each subject. Each sliding window 
% has width winP(1) (default = 30) and slides winP(2) (default = 1) time 
% points at each step.

% vecDFCzcorr is a nPairs x nWin x nSub array with (ijk)th entry 
% corresponding to the Fisher-transformed correlation between region pair i 
% in window j for subject k.

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
if nargin < 2
    winP = [30 1];
end

winLength = winP(1);
winStep = winP(2);
nSub = size(TC,3);
nRegions = size(TC,1);
nT = size(TC,2);

% Get number of windows
win = 1:winStep:(nT+1-winLength);
nWin = length(win);

%% Calculate Dynamic FC
% Set up tapered window
taper_weights = tukeywin(winLength+2);
taper_weights = taper_weights(2:end-1);

% Get indices of upper triangular part of matrix
mat=ones(nRegions,nRegions);
ind = triu(mat,1)==1;

% Preallocate the output
nPairs = nRegions*(nRegions-1)/2;
vecDFCzcorr = zeros(nPairs, nWin, nSub);

for n = 1:nSub
    % Iterate accross windows
    for w1 = 1:nWin
        % Get the current window
        w = win(w1);
        % Get the current data
        winData = squeeze(TC(:,w:w+winLength-1,n));
        
        % Compute the correlation
        DFCcorr = weightedcorrs(winData',taper_weights');
        
        % Fisher transform the correlation matrix
        zCorr = 0.5.* log((1+DFCcorr)./(1-DFCcorr));
        % Vectorise upper triangular part of the correlation matrix
        vecDFCzcorr(:,w1,n) = zCorr(ind);
    end
end
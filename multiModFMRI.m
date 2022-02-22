function [multiModStruct] = multiModFMRI(DFCzcorr, Gamma, Omega)
% Dependency: GenLouvain package - freely available on 
% http://netwiki.amath.unc.edu/GenLouvain/GenLouvain

% Partitions regions into modules using a Louvain-like algorithm to
% maximise multilayer modularity

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% parameter initialisation and preallocation of output
nSub = size(DFCzcorr,3);
nWin = size(DFCzcorr,2);
nPairs = size(DFCzcorr,1);
nRegions = (1+sqrt(1+(8*nPairs)))/2;

% indices of upper triangular part of matrix
mat=ones(nRegions,nRegions);
indTriu = triu(mat,1)==1;

multiModStruct = zeros(nRegions, nWin, nSub);

%% multilayer modularity matrix
for n = 1:nSub
    B = spalloc(nRegions*nWin,nRegions*nWin,nRegions*nRegions*nWin+2*nRegions*nWin);
    twomu = 0;
    for w = 1:nWin
        % recover correlation matrix
        thisCorr = squeeze(DFCzcorr(:,w,n));
        A = zeros(nRegions,nRegions);
        A(indTriu) = thisCorr;
        A = A + A';
        
        % strength of nodes
        k = sum(A);
        twom = sum(k);
        twomu = twomu + twom;
        index = (1:nRegions)+(w-1)*nRegions;
        
        B(index,index) = A - Gamma*(k')*k/twom;
    end

twomu = twomu + 2*Omega*nRegions*(nWin-1);
B = B + Omega*spdiags(ones(nRegions*nWin,2),[-nRegions,nRegions],nRegions*nWin,nRegions*nWin);

%% Louvain-like algorithm
[S,Q] = genlouvain(B,[],0);
Q = Q/twomu;
S = reshape(S,nRegions,nWin);
multiModStruct(:,:,n) = S;

end
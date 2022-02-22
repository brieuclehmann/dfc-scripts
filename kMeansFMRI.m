function [IDXall, centroids, sum] = kMeansFMRI(data, nClust)
% kMeansFMRI runs k-means clustering algorithm on pooled data for 
% k = 1,...,nClust. IDXall contains the centroid label 

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
% Set default number of clusters
if nargin < 2
    nClust = 12;
end

nSub = size(data,3);
nWin = size(data,2);
nPairs = size(data,1);

IDXall = zeros(nClust, nSub, nWin);
centroids = zeros(nClust, nClust, nPairs);
sum = zeros(1,nClust);

%% Run k-means clustering
for k = 1:nClust
    k
    % Apply k-means algorithm
    [thisIDX, C, thissum] = kmeans(data(:,:)',k,'distance', 'sqeuclidean', 'Replicates', 40);
    % Store centroids
    centroids(k,1:k,:) = C;
    % Store cluster indices for each subject in each window
    IDXall(k,:,:) = reshape(thisIDX, size(data,2), size(data,3))';
    % Store 2-norm of sums of point-to-centroid distances (for elbow plot)
    sum(k) = norm(thissum);
end
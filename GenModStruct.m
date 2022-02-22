function modStruct = GenModStruct(nMod, nRegions, nStates, nmiVal)
% Generates a module structure for nStates, partitioning nRegions ROIs into
% exactly nMod modules. nmiVal determines the maximum normalised mutual
% information allowed between any two states (i.e. it controls how similar 
% any two states can be). modStruct is a nStates x nRegions with (ij)th 
% entry corresponding to the module assignment of region j in state i

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
% Set default nmi threshold for maximum similarity between states
if nargin < 4
    nmiVal = 0.5;
end

maxnmi = 1;

%% Generate module structure
while maxnmi >= nmiVal % Check nmi threshold
    modStruct = zeros(nStates,nRegions);
    % Assign a module structure to each state
    for state = 1:nStates
        while length(unique(modStruct(state,:))) < nMod
            % Ensure there are exactly nMod modules in each state
            modStruct(state,:) = randi(nMod,1,nRegions);
        end
    end
    
    % Check similarity of states by calculating nmi between each pair
    nmiStates = zeros(nStates);
    for s1=1:nStates
        for s2=1:nStates
            nmiStates(s1,s2)=nmi(modStruct(s1,:), modStruct(s2,:));
        end
    end
    % Set diagonal to zero to ignore nmi of a state with itself
    nmiStates(eye(size(nmiStates))~=0) = 0;
    % Find value of nmi for most similar pair
    maxnmi = max(max(nmiStates));
end
function stateSeq = GenStateSeq(nStates, varargin)
% GenStateSeq generates state sequences of length nDiffStates (default = 4)
% from nStates attainable states for nSub (default = 1) subjects. For each
% subject, the state sequence is generated as a nStates-state Markov chain
% with transition matrix transition_Mat. The default transition_Mat assigns
% an equal probability to transition to another state, and does not allow
% transitions back to the same state. 

% stateSeq is a nDiffStates x nSub matrix with (ij)th entry corresponding
% to the (i)th state of subject j. 

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
parseObj = inputParser;
% Default number of states traversed (i.e. number of state changes + 1)
default_nDiffStates = 4;
addOptional(parseObj,'nDiffStates',default_nDiffStates,@isnumeric);
% Default number of subjects in group
default_nSub = 1;
addOptional(parseObj,'nSub',default_nSub,@isnumeric);
% Default state-transition matrix
default_transitionMat = (1/(nStates-1)).*ones(nStates,nStates);
default_transitionMat(eye(nStates)~=0) = 0;
addOptional(parseObj,'transitionMat',default_transitionMat,@ismatrix);
parse(parseObj, varargin{:});

nDiffStates = parseObj.Results.nDiffStates;
nSub = parseObj.Results.nSub;
transitionMat = parseObj.Results.transitionMat;

% Calculate equilibrium distribution
[V, D] = eig(transitionMat','vector');
[~, ind] = sort(D,'descend');
initDist = V(:,ind)/sum(V(:,ind));


%% generate a state sequence
stateSeq = zeros(nDiffStates, nSub);
% Pick an initial state based on equilibrium distribution
r = rand(1,nSub);
stateSeq(1,:) = arrayfun(@(z)sum(z >= cumsum([0, initDist'])), r);

for n = 1:nSub
    for i = 1:nDiffStates-1
        r = rand;
        s = 1;
        p = transitionMat(stateSeq(i,n),s);
        % Choose the next state based on the previous state
        while r > p
            s = s + 1;
            p = p + transitionMat(stateSeq(i,n),s);
        end
        stateSeq(i+1,n) = s;
    end
end
function [eventSeq, moduleEvents] = GenEventSeq(stateSeq, modStruct, varargin)
% Generates a sequence of events of length nT for nRegions ROIs comprised
% of a module-specific event sequence and a region-specific event
% sequence. A module-specific event sequence is a two-state Markov chain 
% with equilibrium probability pMod (default = 0.5) and lag-1 
% autocorrelation moduleAC (default = 0). A region-specific event sequence
% is a two-state Markov chain with equilibrium probability pReg 
% (default = 0.5) and lag-1 autocorrelation regionAC (default = 0). 
% Module-specific events have an amplitude of aModule (default = 2) 
% relative to region-specific events. 

% eventSeq is a nRegion x nT x nSub array with (ijk)th entry corresponding 
% to the number of events in region i for subject k at time point j

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
nSub = size(stateSeq,2); 
nDiffStates = size(stateSeq,1);
nMod = max(max(modStruct));
nRegions = size(modStruct,2);

parseObj = inputParser;
% Default amplitude of module-specific event relative to region-specific 
% event 
default_aModule = 2;
addOptional(parseObj,'aModule',default_aModule,@isnumeric);
% Default equilibrium probability of region-specific events
default_pRegionEvent = 0.5;
addOptional(parseObj,'pRegionEvent',default_pRegionEvent,@isnumeric);
% Default autocorrelation of region-specific events
default_regionAC = 0;
addOptional(parseObj, 'regionAC',default_regionAC,@isnumeric);
% Default equilibrium probability of module-specific events
default_pModuleEvent = 0.5;
addOptional(parseObj,'pModuleEvent',default_pModuleEvent,@isnumeric);
% Defaut autocorrelation of module-specific events
default_moduleAC = 0;
addOptional(parseObj, 'moduleAC',default_moduleAC,@isnumeric);
% Default number of time points
default_nT = 360;
addOptional(parseObj,'nT',default_nT,@isnumeric);
parse(parseObj, varargin{:});

aModule = parseObj.Results.aModule;
pRegionEvent = parseObj.Results.pRegionEvent;
pModuleEvent = parseObj.Results.pModuleEvent;
regionAC = parseObj.Results.regionAC;
moduleAC = parseObj.Results.moduleAC;
nT = parseObj.Results.nT;

% Region-specifc event probabilities
pReg = [regionAC + pRegionEvent*(1-regionAC), pRegionEvent*(1-regionAC)];

% Module-specific event probabilities
pMod = [moduleAC + pModuleEvent*(1-moduleAC), pModuleEvent*(1-moduleAC)];

% Check parameters yield valid Markov chain
if any(pMod <= 0) || any(pMod >= 1)
    warning('Must have -moduleAC < pModuleEvent*(1-moduleAC) < 1')
end

if any(pReg <= 0) || any(pReg >= 1)
    warning('Must have -regionAC < pRegionEvent*(1-regionAC) < 1')
end

%% Time spent in each state
tState = floor(nT/nDiffStates).*ones(1,nDiffStates);
i = 1;
% Ensure that time between state changes is as constant as possible
while sum(tState) < nT
    tState(i) = tState(i) + 1;
    i = i + 1;
end

%% Event generation
% Generate module-specific events

moduleEvents = zeros(nRegions, nT, nSub);
regionEvents = zeros(nRegions, nT, nSub);

for n = 1:nSub
    % Get this subject's state sequence
    thisSeq = repelem(stateSeq(:,n),tState);
    
    % Initial events
    % Generate intial module events
    modEvents = binornd(1,pModuleEvent,1,nMod);
    
    % Get the inital module structure
    thisMod = modStruct(thisSeq(1),:);
    % Assign initial module-specific events for each region
    moduleEvents(:,1,n) = modEvents(thisMod);
    % Generate initial region-specific events
    regionEvents(:,1,n) = binornd(1,pRegionEvent,1,nRegions);
    
    for i = 2:nT
        % Generate module events at time i
        for mod = 1:nMod
            if rand < pMod(2-modEvents(mod))
                modEvents(mod) = 1;
            else
                modEvents(mod) = 0;
            end
        end
        
        % Get the current state
        thisState = thisSeq(i);
        % Get module structure in this state
        thisMod = modStruct(thisState,:);
        % Assign module-specific events for each region
        moduleEvents(:,i,n) = modEvents(thisMod);
        
        % Generate region-specific events at time i
        for reg = 1:nRegions
            if rand < pReg(2-(regionEvents(reg,i-1,n)~=0));
                regionEvents(reg,i,n) = 1;
            else
                regionEvents(reg,i,n) = 0;
            end
            
        end
    end
end

% Add region-specific events array to amplitude-adjusted
% module-specific events array to yield total sequence of events
eventSeq = aModule*moduleEvents + regionEvents;
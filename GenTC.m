function TC = GenTC(eventSeq, varargin)
% Dependency: SPM
% GenTC generates time series for each subject by convolving a region's 
% event sequence eventSeq with a haemodynamic response function 
% (default = SPM default) at a repetition time of TR (default = 2s), 
% rescaling to have standard deviation 1, adding white noise of standard 
% deviation aNoise (default = 0.2), and applying a high-pass filter.

% TC is a a nRegion x nT x nSub array with (ijk)th entry corresponding 
% to the simulated BOLD response of region i for subject k at time point j

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
parseObj = inputParser;
% Default SD of white noise (implies default SNR = 5)
default_aNoise = 0.2;
addOptional(parseObj,'aNoise',default_aNoise,@isnumeric)
% Default repetition time
default_TR = 2;
addOptional(parseObj,'TR',default_TR,@isnumeric);
% Default HRF parameters - (use SPM default)
default_P = [];
addOptional(parseObj,'P',default_P,@isvector);
% Apply high-pass filtering? (default = yes)
default_filter = 1;
addOptional(parseObj,'filter',default_filter,@islogical);
parse(parseObj, varargin{:});

TR = parseObj.Results.TR;
P = parseObj.Results.P;
aNoise = parseObj.Results.aNoise;
filter = parseObj.Results.filter;

nSub = size(eventSeq,3);
nRegions = size(eventSeq,1);
nT = size(eventSeq,2);

% Get HRF function
HRF = spm_hrf(TR,P);

% High-pass filter parameters
param.dt = 2;
param.filt_fg = 1/30; % for window length 30
ny = 1/2/param.dt;
Wn = param.filt_fg./ny;
[filtb,filta] = butter(4,Wn,'high');

%% Generate time courses
TC = zeros(nRegions, nT, nSub);

for n = 1:nSub
    for i = 1:nRegions
        % Convolve the event sequence with the HRF
        thisTC = conv(HRF,eventSeq(i,:,n));
        TC(i,:,n) = thisTC(1:nT);
        % Rescale to have SD equal to 1
        TC(i,:,n) = TC(i,:,n)/std(TC(i,:,n));
        % Add some random noise (SNR = 1/aNoise)
        TC(i,:,n) = TC(i,:,n) + aNoise*randn(size(TC(i,:,n)));
        
        if filter
            % Apply high pass filter to remove any frequencies below window size
            TC(i,:,n) = filtfilt(filtb,filta,TC(i,:,n));
        end
    end
end
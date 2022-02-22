function [err_struct, err_flex, mmFlex, trueFlex] = errorMultiMod(multiModStruct, winModStruct)
% Computes error measures for multilayer modularity analysis. err_struct is 
% the percentage connectivity, err_flex is the error in flexibility, mmFlex 
% is the recovered flexibility, and trueFlex is the true flexibility.

%   Author: Brieuc Lehmann
%   E-mail: brieuc.lehmann@mrc-bsu.cam.ac.uk
%   Date: 3 April 2017

%% Parameter initialisation
nSub = size(multiModStruct,3);
nWin = size(multiModStruct,2);
nRegions = size(multiModStruct,1);

err_struct = zeros(nSub,nWin);
err_flex = zeros(nSub,1);
mmFlex = zeros(nSub,1);
trueFlex = zeros(nSub,1);

%% Computation of error measures
for n = 1:nSub
    % Compute flexibility measures
    trueFlex(n) = sum(sum(diff(winModStruct(:,:,n),1,2) ~= 0))/(nRegions*(nWin - 1));
    mmFlex(n) = sum(sum(diff(multiModStruct(:,:,n),1,2) ~= 0))/(nRegions*(nWin - 1));
    err_flex(n) = (mmFlex(n) - trueFlex(n));
    
    for w = 1:nWin
        % Get recovered module structure
        thisStruct = squeeze(multiModStruct(:,w,n));
        % Get true module structure
        trueStruct = squeeze(winModStruct(:,w,n));
       
        % Compute percentage error in connectivity structure
        thisConnMatrix = repmat(thisStruct,1,length(thisStruct)).^2 == thisStruct*thisStruct';
        trueConnMatrix = repmat(trueStruct,1,length(trueStruct)).^2 == trueStruct*trueStruct';
        err_struct(n,w) = sum(sum(thisConnMatrix ~= trueConnMatrix))/(nRegions^2);
    end
end
% Take mean (across windows) error in connectivity structure
err_struct = mean(err_struct,2);
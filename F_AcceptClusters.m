function [gm2] = F_AcceptClusters(gm, nsubcl, majorN)
global CFG
% F_AcceptClusters
% To accept clusters that are representative for the majority of participants.
%
% INPUT
% gm -  gaussian mixture from F_GM2 cell 1 x nroi. Variable number of
%       gaussian components per ROI.
% nsubcl - number of subjects that contributed to particular gaussian component.
% majorN - number of subjects that have to contribute to particular
%           component to be accepted.
%
% OUTPUT
% gm2 - gaussian mixture from F_GM2 cell 1 x nroi. Variable number of
%       gaussian components per ROI. 
%
%
gm2 = cell(1, CFG.nroi);
for iroi = 1:CFG.nroi
    if isempty(gm{iroi})
        continue;
    end
    strongclustidx = find(nsubcl{iroi} >= majorN);
% too much inputs in gmdistribution constructor
%     gm2{iroi} = gmdistribution( 'ComponentProportion', gm{iroi}.ComponentProportion , ...	
%                                 'CovarianceType', gm{iroi}.CovarianceType, ...
%                                 'DistributionName'	,gm{iroi}.DistributionName, ...	
%                                 'mu', gm{iroi}.mu(strongclustidx,:), ...	
%                                 'NumComponents'	, length(strongclustidx), ...	
%                                 'NumVariables'	, gm{iroi}.NumVariables, ...	
%                                 'ProbabilityTolerance'	, gm{iroi}.ProbabilityTolerance, ...	
%                                 'SharedCovariance', gm{iroi}.SharedCovariance, ...	
%                                 'Sigma', gm{iroi}.Sigma(:, :, strongclustidx));
    gm2{iroi} = gmdistribution( gm{iroi}.mu(strongclustidx,:), ...	 
                                gm{iroi}.Sigma(:, :, strongclustidx), ...
                                gm{iroi}.ComponentProportion(strongclustidx));
end
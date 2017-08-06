function [allrois, subindex] = F_Rearrange(mode)
global CFG
% Loads and gathers 1stGM data from individuals to one structure (this is not  pooling, but just rearranging structures)
% INPUT
% mode -  determines file with 1GM data to load 'main' - GM models for generating spectral profiles, 'class' - temporary GM models for do classification iterative procedure (Fig 3a);
%
% OUTPUT
% allrois - cell 116x1 nclust x nfreq x nsub
% subindex - [1 1 1 ... 1 2 2 2 ... 2 3 3 3 ... 3 ...22 22 22 ... 22]'

% !!! It is very important not to mix the order of 1st-lvl clusters when
% rearranging - otherwise the time duration information will be disturbed
    allsubs = zeros(CFG.nclust, CFG.nfreq, CFG.nroi, CFG.nsub);
    for isub = 1:CFG.nsub 
        if(strcmp(mode,'main'))
            load(sprintf('%s%s_nclust%d_sub%d', CFG.rsltsDir, CFG.indivname, CFG.nclust, CFG.goodsub(isub)));  % roicluster - nclust x nfreq x nroi
            allsubs(:,:,:,isub) = roicluster;
        elseif(strcmp(mode,'class'))
            load(sprintf('%s%s_testdata%d_sub%d', CFG.rsltsDir, CFG.indivname, CFG.nclust, CFG.goodsub(isub)));  % roicluster - nclust x nfreq x nroi
            allsubs(:,:,:,isub) = roicluster;
        end
    end
    allrois = cell(CFG.nroi,1); % 116x1 cell array: 10 clusters x 42 freqs x 22 subjects double    


    for iroi = 1:CFG.nroi,
        allrois{iroi}(:,:,:) = squeeze(allsubs(:,:,iroi,:));
    end
    clear allsubs

    % build sub index to keep track of subs
    subindex = [];
    for isub = 1:CFG.nsub
        subindex = [subindex ones(1, CFG.nclust)*CFG.goodsub(isub)];
    end
    subindex = subindex';

end
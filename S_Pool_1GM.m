% S_Pool_1GM

% load first-level clusters and store in one matrix
allsubs = zeros(CFG.nclust, CFG.nfreq, CFG.nroi, CFG.nsub);
for isub = CFG.goodsub
    %filename = sprintf('%s%s_nclust%d_sub%d',HomeFolder,indivname,nclust,isub);
    filename = sprintf('%s%s_nclust%d_sub%d', CFG.rsltsDir, CFG.indivname, CFG.nclust, isub);
    load(filename,'roicluster','ff');
    allsubs(:,:,:,isub) = roicluster;
end
% optional: remove all empty entries, i.e. if you have missing subjects
% e.g., allsubs(:,:,:,badsub)=[];

% rearrange data
allrois = cell(nroi,1);
for iroi = 1:nroi,
    allrois{iroi}(:,:,:) = squeeze(allsubs(:,:,iroi,:));
end
clear allsubs

% build sub index to keep track of subs
subindex = [];
for isub = 1:CFG.nsub,
    subindex = [subindex ones(1, CFG.nclust)*isub];
end
subindex = subindex';
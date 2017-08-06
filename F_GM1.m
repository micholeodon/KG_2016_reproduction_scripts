function [roicluster, pt] = F_GM1(ROIidx, sourcef)
global CFG
% F_GM1
%

% ROIidx - list of ROI's indices to analyse
%
% OUTPUT:
% roicluster - matrix nclust x nfreq x nroi with GM centroids in R^{nfreq} space
% pt - cluster duration (number of trials in K-means cluster = number of (not subsequent!) seconds)

    disp('STARTING 1st lvl GM modelling ...')
    disp('Loading sourcemodel ...')
    load(CFG.SourceMod);    % [sourcemodel]
    roicluster = zeros(CFG.nclust, CFG.nfreq, CFG.nroi);          % create matrix for individual cluster
    id = cell(CFG.nroi, 1);                              % create matrix for cluster indices
    pt = zeros(CFG.nroi, CFG.nclust);                        % create matrix for cluster durations
    opts = statset('Display','final');
    % CLUSTERING whos
    % average each atlas ROI and cluster power spectra (kmeans and GM models)
    disp('ROI loop ...')
    for iroi = ROIidx
        iroi
        % find all voxels belonging to a ROI
        roivoxel = find(CFG.sourceatlas.tissue==iroi);

        % match the voxel indices between the sourcemodel (which includes inside
        % and outside voxels) with our voxel indices (which only has inside voxels)
        [~,insideroivoxel] = intersect(sourcemodel.inside,roivoxel);

        roimean = squeeze(mean(sourcef(:,insideroivoxel,:),2)); % sourcef is global
        disp('removing outliers ...')
        %------------------------------
        olier = zscore(mean(roimean,2));            % clean data in frequency domain
        fi = find(olier>2.5);                       % by removing outlier
        roimean(fi,:) = [];                         % (optional)
        %------------------------------

        roimean = roimean-1;                    % make values below/above zero (instead of one)

        % k-means clustering
        disp('K-means ...')
        [idd, rr] = kmeans(roimean, CFG.nclust,'Options', opts,'Distance','cosine','Replicates', 10);
        
        disp(iroi) % TEST !!!
        disp(size(roicluster)) % TEST !!!
        % Gaussian Mixture clustering; uses k-means cluster as starting point
        disp('Gaussian modelling ...')
        gm = gmdistribution.fit(roimean, CFG.nclust, 'Start', idd, 'Options', opts, 'Regularize', CFG.regu);
        roicluster(:,:,iroi) = gm.mu;           % 1-st level cluster
        id{iroi} = cluster(gm,roimean);         % index cluster

        disp('Computing trial duration ...')
        %-----------------------------------------------
        for iclust = 1:CFG.nclust                                   % get the amount of trials belonging
            pt(iroi,iclust) = length(find(id{iroi}==iclust));   % to each cluster (expressed in %)
            pt(iroi,iclust) = pt(iroi,iclust)/length(id{iroi})*100;
        end
        %-----------------------------------------------

        
    end



end
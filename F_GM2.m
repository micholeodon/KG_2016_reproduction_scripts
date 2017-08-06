function [gm, nsubcl, allidx] = F_GM2(allrois, optopt, numclust, subjIdx, subindex)

    global CFG
    % INPUT
    % allrois - cell 1x116 (nclust, nfreq, nsubj) rearranged data
    % optopt - ['fix' | 'opt'] - string; decides about what form numclust
    % input variable should be
    % numclust - if optopt = 'fix': 1x1 integer > 0; if optopt = 'opt' it loads optimal cluster information from optclust.mat file CFG.optname found in CFG.rsltsDir
    % subjIdx - vector of subject indices to be pooled into 2nd lvl GM model
    % subindex - variable from F_Rearrange() to keep track of subs

    % OUTPUT
    % gm - cell 116x1 (gaussian model for each cluster with different or fixed numbr of centroids per roi)
    % nsubcl -  cell 116x1 (nclust) subject contribution for each cluster (number of subjects, not their indices !!!). 
    %           Length of each vector is equal to the number of optimal clusters ! 
    % allidx -  indices of cluster in 2nd-level models (not 1st because, there they consist numbers from 0-15, 
    %           so they correspond to optimal number of clusters range, not to the fixed 10 clusters in 1st GM)
    %           for each point for each roi

disp('Entering goodroi loop ...')
    for iroi = CFG.goodroi
        iroi
        % POOLING - rearrange data (all 11 subject clusters in 1st dimension, frequencies in 2nd dimension)
        disp('Pooling, extracting and rearranging data ...') % extracting = taking subset of allrois corresponding with indices in subjIdx
        ss = find(ismember(CFG.goodsub, subjIdx)); % which subject fibre extract from allrois{roi}
        roidata = allrois{iroi}(:,:,ss);        
        roidata = permute(roidata(:,:,:),[1 3 2]); % trick to do reshape correctly
        roidata = reshape(roidata, CFG.nclust*length(subjIdx), CFG.nfreq); % if 22 total good subjects and 10 clusters in 1stGM then roidata should be 110x42
        
        disp('Loading info about optimal number of clusters ...')
        if(strcmp(optopt,'fix'))
            optn = numclust;
        elseif(strcmp(optopt,'opt'))
            load(CFG.optname,'nclustopt')
            optn = nclustopt(iroi)
        end
        
        % k-means with optimal number of clusters
        disp('k-means ...')
        opts = statset('MaxIter',1000);             % maximum number of iterations
        [idd,rr] = kmeans(roidata, optn,'Options', opts,'Distance','cosine','Replicates',10);
        
        % Gaussian Mixture clustering; uses k-means cluster as starting point
        disp('gaussian modelling ...')
        opts = statset('Display','final');
        gm{iroi} = gmdistribution.fit(roidata, optn, 'Start', idd, 'Options', opts, 'Regularize', CFG.regu);
        idx = cluster(gm{iroi}, roidata); % cluster indices for each of pooled points
        allidx(:,iroi) = idx;               % indices of cluster in 1st-level models for each roi
        
        %-----------------------------------
        disp('Calculating subject contribution ...')
        for k = 1:optn,                     % number of subjects contributing to each of 2nd lvl clusters
            fi = find(idx==k);              % (important if you only want to show cluster of majority of participants)
            nsubcl{iroi}(k) = length(unique(subindex(fi)));
        end
        %-----------------------------------
        
        disp([num2str(iroi),' DONE !']);
    end


end
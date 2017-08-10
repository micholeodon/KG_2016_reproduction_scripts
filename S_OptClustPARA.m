% Evaluate optimal number of clusters in 2nd lvl GM Model

%optname = sprintf('%soptclust',HomeFolder);
optname = CFG.optname


% MICHAK : to run optimal number of clusters calculation in parallel
% Run 4 Matlab shells in the same time (so optimal number of clusters will be divided to four files. Don't worry it will be merge by S_SAVEOptClustPARA script.)
% Your job here is to change multiplicant to 0, 1, 2, and finally 3 each time you run this. Same in niterF but 1,2,3 and 4. Nothing more to change.
% After you run this script 4 times with changing this numbers, then you should run S_SAVEOptClustPARA to merge results and save.
%  0, 1, 2, 3
iiter0 = 1*(CFG.niter/4) + 1;
%  1, 2, 3, 4
niterF = 2*(CFG.niter/4);

lettset = {'A','B','C','D'};  
%   'A', 'B', 'C',  'D'
lett = lettset{4*niterF/CFG.niter}

allnclust.(lett) = NaN(CFG.nroi, CFG.niter);    % create matrix for optimal cluster numbers
allsilhval.(lett) = NaN(CFG.nroi, CFG.niter);   % create matrix for Silhoutte values

for iiter = iiter0:niterF
    for iroi = CFG.goodroi
        iiter % michak: added for progress tracking
        iroi % michak: added for progress tracking
             % rearrange data
        roidata = allrois{iroi};
        roidata = permute(roidata,[1 3 2]);
        roidata = reshape(roidata, CFG.nclust*CFG.nsub, CFG.nfreq);
        tmp = roidata;
        
        % evaluate cluster solutions
        E = evalclusters(tmp,'kmeans','silhouette','klist', CFG.clustnum,'Distance','cosine');
        allnclust.(lett)(iroi,iiter) = E.OptimalK;                     % store optimal number for each iteration
        allsilhval.(lett)(iroi,iiter) = E.CriterionValues(E.OptimalK); % store Silhouette value for each iteration
    end
    save([CFG.rsltsDir, 'tmp_allclustsilh', '_', num2str(iiter0), '_', num2str(niterF)], 'allnclust', 'allsilhval', 'iiter', 'iroi', 'iiter0', 'niterF') % michak: added in case of need of stopping and not losing partial results
end

%

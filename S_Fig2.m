%% ROI spectral profiles (spectrum, surf)
% We've used a modified version of the mseb.m function
% (multiple shaded error bars) to plot modes and SEM
% http://www.mathworks.com/matlabcentral/fileexchange/47950-mseb-x-y-errbar-lineprops-transparent-



mkdir([CFG.rsltsDir, 'fig/'])

roi = 5; % ROI you want to plot, 9 is Precentral Gyrus Left; check for example 1, 5, 17.

load([CFG.rsltsDir, 'ff_all'], 'ff_all') % for later plotting ff for individual
load(CFG.finalname, 'gm', 'nsubcl', 'allidx', 'subindex', 'gm_ok')

% select clusters representative for majority of participants
disp('Selecting cluster present in majority of participants ...')
strongclustidx  = sort(find(nsubcl{roi} >= CFG.majorN))
numsubjpercomp  = [nsubcl{roi}(strongclustidx)];
clusters        = gm{roi}.mu; % here are already cluster truncated with F_AcceptClusters.m

disp('Estimating cluster duration ...')
whos allidx
cldur = nan(1, length(strongclustidx)); % preallocate
for c2 = strongclustidx 
    c2 % cluster index
    pidx = find(ismember(allidx(:,roi), c2))' % indices of points belonging to particular 2nd lvl cluster c2
    subjlist = unique(subindex(pidx))' % get subject indices to know which pt should be loaded (pt tells how long each of 1st-lvl centroid lasts during recording)
    
    all_pt_roi = nan(CFG.nsub, CFG.nclust); % preallocate; matrix subject x 1-stlvl cluster which shows how much each cluster lasts in current roi
   
    % get durations for each subject and cluster
    for ss = subjlist
        load([CFG.rsltsDir, CFG.indivname, '_nclust10_sub', num2str(ss)], 'pt')
        all_pt_roi(ss,:) = pt(roi,:)
    end
    
    % estimate cluster duration
    pdur = nan(1, CFG.nsub*CFG.nclust);
    for pp = pidx
        
        % get 1st-lvl cluster index
        c1 = mod(pp, CFG.nclust); % get last digit 
        if c1 == 0
            c1 = 10; % zero means cluster no 10
        end
        c1
        
        ss = subindex(pp)
        pdur(pp) = all_pt_roi(ss, c1)
    end
    
    % Decide which is correct duration time estimation:
    % 1. Get mean duration of all points contributing to cluster
    %cldur(find(ismember(strongclustidx, c2))) = nanmen(pdur) % nenmean ? max ? 
    % 2. Get max duration of all points contributing to cluster
    %cldur(find(ismember(strongclustidx, c2))) = max(pdur) 
    % 3. Sum all points durations (= sum number of trials that each cluster was evident) and divide by the number of subject that
    % contributed to this cluster ("averaging across subject" mentioned in
    % the article). This method gives percentages very similar to those
    % from article !!!
    cldur(find(ismember(strongclustidx, c2))) = nansum(pdur)/nsubcl{roi}(c2) 
end



disp('Estimating cluster variance ...')
Ncl             = gm{roi}.NumComponents;
clusters_std    = zeros(Ncl,CFG.nfreq); % for each of 42 freqs we should calculate SE across subjects (need to go back for 1GM model)
for cl = 1:Ncl
    %     clusters_std(cl,:) = std(clusters(cl,:))/sqrt(22); % WRONG !!!!
    % Instead I will take diagonal elements of covariance matrix and divide by sqrt of
    % number of subjects that contributed to this. This is maybe not exact but
    % but more rigorous than division by sqrt of number of points that
    % contributed for a component (ten times smaller number, because each
    % subject gave 10 centroids).
    cl
    var_est             = diag(gm{roi}.Sigma(:,:,cl));
     
    Nsubj               = numsubjpercomp(cl);
    SE                  = sqrt(var_est)./sqrt(Nsubj);
    clusters_std(cl,:)  = SE;
end
disp('Variance estimated !')

% Generate figure
figure, h = mseb(ff_all{1}, clusters, clusters_std);
title(sprintf('Raw spectral modes in ROI %s: %s', num2str(roi), CFG.atlas.tissuelabel{roi}),'Interpreter','none','FontSize',12,'FontWeight','bold');
xlabel('Frequencies (Hz)','FontSize',12);
ylabel('Normalised power','FontSize',12);
set(gca,'XGrid','on','XMinorTick','on','XScale','log','XTick',[1 3 5 10 20 50 81],'XTickLabel',{'1','3','5','10','20','50','81'});
lgd = legend(num2str(round(cldur)'));
lgd.String = strcat(lgd.String, ' %');
xlim([1 120])

% saving figure
saveas(gcf,[CFG.rsltsDir, 'fig/roi', num2str(roi), '.png'])


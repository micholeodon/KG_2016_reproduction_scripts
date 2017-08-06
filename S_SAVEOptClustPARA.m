% AFTER ALL: just load and merge variables into 
% save(optname,'allnclust','allsilhval');
% Save optimal no of clusters for each iteration plus silhouette values.
% Cosine distance criterium is important (others don't capture shape of
% amplitudes so well).


JJ1 = [0*(CFG.niter/4) + 1 , 1*(CFG.niter/4)] % eg. for niter=1000  [1 250]
JJ2 = [1*(CFG.niter/4) + 1 , 2*(CFG.niter/4)] % eg. for niter=1000  [251 500]
JJ3 = [2*(CFG.niter/4) + 1 , 3*(CFG.niter/4)] % eg. for niter=1000  [501 750]
JJ4 = [3*(CFG.niter/4) + 1 , 4*(CFG.niter/4)] % eg. for niter=1000  [751 1000]

clA = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ1(1)) '_' num2str(JJ1(2)) '.mat'], 'allnclust')
clB = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ2(1)) '_' num2str(JJ2(2)) '.mat'], 'allnclust')
clC = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ3(1)) '_' num2str(JJ3(2)) '.mat'], 'allnclust')
clD = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ4(1)) '_' num2str(JJ4(2)) '.mat'], 'allnclust')
allnclust = nan(CFG.nroi, CFG.niter); 
allnclust(:, JJ1(1):JJ1(2)) = clA.allnclust.A(:,JJ1(1):JJ1(2));
allnclust(:, JJ2(1):JJ2(2)) = clB.allnclust.B(:,JJ2(1):JJ2(2));
allnclust(:, JJ3(1):JJ3(2)) = clC.allnclust.C(:,JJ3(1):JJ3(2));
allnclust(:, JJ4(1):JJ4(2)) = clD.allnclust.D(:,JJ4(1):JJ4(2));
siA = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ1(1)) '_' num2str(JJ1(2)) '.mat'], 'allsilhval')
siB = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ2(1)) '_' num2str(JJ2(2)) '.mat'], 'allsilhval')
siC = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ3(1)) '_' num2str(JJ3(2)) '.mat'], 'allsilhval')
siD = load([CFG.rsltsDir,'tmp_allclustsilh_' num2str(JJ4(1)) '_' num2str(JJ4(2)) '.mat'], 'allsilhval')
allsilhval = nan(CFG.nroi, CFG.niter); 
allsilhval(:, JJ1(1):JJ1(2)) = siA.allsilhval.A(:,JJ1(1):JJ1(2));
allsilhval(:, JJ2(1):JJ2(2)) = siB.allsilhval.B(:,JJ2(1):JJ2(2));
allsilhval(:, JJ3(1):JJ3(2)) = siC.allsilhval.C(:,JJ3(1):JJ3(2));
allsilhval(:, JJ4(1):JJ4(2)) = siD.allsilhval.D(:,JJ4(1):JJ4(2));

% get optimal cluster numbers from Silhoutte values
nclustopt = NaN(CFG.nroi,1);                    % create matrix for optimal cluster numbers
[~,iclustmaxid] = max(allsilhval,[],2);
for iroi = CFG.goodroi,
    nclustopt(iroi) = allnclust(iroi, iclustmaxid(iroi));
end;


optname = CFG.optname;
save(optname, 'allnclust','allsilhval','nclustopt');

clear JJ1 JJ2 JJ3 JJ4
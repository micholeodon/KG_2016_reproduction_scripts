Z = linkage(roifit(CFG.goodroi, CFG.goodroi), 'average', 'cosine');
dendrogram(Z,0)
[~, T] = dendrogram(Z,20); % tree displayed starting from 20 clusters
size(T)

% generate table with tissue labels
rows = [];
for tt = 1:20
  rows.idx{tt} = find(T==tt);
  rows.lab{tt} = {CFG.atlas.tissuelabel{rows.idx{tt}}}';
end
rows
rows.lab
rows.idx
rows.idx{:}
rows.lab{:}

%

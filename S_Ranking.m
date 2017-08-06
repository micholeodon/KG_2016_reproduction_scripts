for jj = 1:CFG.Niter
    ranks(:,jj) = diag(iterdata(jj).roifitrank); % diag because we consider recognizing correct ROI (match)
end
% - trim 20% lowest values
% - trim 20% highest values
% - compute mean rank
% ... all in one line :)
t20meanrank = trimmean(reshape(ranks,1,[]),20); % one general rank; trimmean ignores NaN values

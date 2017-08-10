% This is only valid if you use AAL atlas, because labels from 1 to 108 are
% orderered pairwaise: odd labels are L-structures, even labels are
% R-structures.
% Labels from 109-116 are Vermis labels, which are central structure
% without division to lateral homologues.
% Check this here and you can confirm location of each structure with S{1..10}_Fig.TIF :
% You can also check reindex.txt file to see correspondence with atlas numeration of structures and order of structures in authors' supporting materials (S1-S10).
{CFG.atlas.tissuelabel{1:2:116};CFG.atlas.tissuelabel{2:2:116}}


for jj = 1:CFG.Niter
    
    % for ROI with homologues use this procedure ...
    for rr = 1:108
        if(mod(rr,2)) 
            % take next
            that = iterdata(jj).roifitrank(rr,rr);
            next = iterdata(jj).roifitrank(rr,rr + 1);
            R = min(that, next); % take lower rank from two homologue ROI
        else
            % take prev
            that = iterdata(jj).roifitrank(rr, rr);
            prev = iterdata(jj).roifitrank(rr, rr - 1);
            R = min(that, prev); % take lower rank from two homologue ROI
        end
        ranks_homo(rr,jj) = R;
    end
    
    % for ROI without homologues use procedure from S_Ranking.m
    ranks_homo(109:116, jj) = diag(iterdata(jj).roifitrank(109:116, 109:116)); % diag because we consider recognizing correct ROI (match)
    
end

% - trim 20% lowest values
% - trim 20% highest values
% - compute mean rank
% ... all in one line :)
t20meanrank_homo = trimmean(reshape(ranks_homo,1,[]),20); % one general rank; trimmean ignores NaN values.

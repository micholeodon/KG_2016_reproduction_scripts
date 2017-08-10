function [roifit, roifitrank] = F_NlogLRank(allrois, gm, subjIdx, mode)
global CFG
% INTPUT
% allrois       - 116x1 cell: nclust x nfreq x nsub
% gm            - gaussian model (2nd lvl) from training set
% subjIdx       - indices of subjects whos 1st GM are to be considered
% mode          - string. 'pool' - compute roifit/roifitrank matrices on pooled data
%                 from all subjects. 'subj_avg' - compute roifit/roifitrank
%                 matrices for each subject separately, then average across
%                 subjects.

% OUTPUT
% roifit - matrix (row-wise) nclust x nclust with negative log-likelihood of fitting data from 1st GM to model 2nd GM
% roifitrank - ranked version (row-wise) of roifit - smaller the rank lower the nlogl better the fit

% !!! Important: roifit and roifit matrices should have meaningful values
% only in entries CFG.goodroi, CFG.goodroi (if some roi are bad they will
% get NaN value).

    switch mode
        case 'pool'
            disp('POOL mode !!!')
            % CALCULATE POSTERIORS
            disp('CALCULATING POSTERIORS and Negative Log-likelihood. NlogL will be ranked additionaly.')

            roifit      = nan(CFG.nroi,CFG.nroi);
            roifitrank  = nan(CFG.nroi,CFG.nroi);

            
            disp('Entering goodroi loop ...')
            iroi = 1;
            jroi = 1;
            for iroi = CFG.goodroi
                iroi

                % pool data
                % POOLING - rearrange data (all 11 subject clusters in 1st dimension, frequencies in 2nd dimension)
                % (Pooling is my assumption, because they have not provided this detail in the paper.)
                disp('Pooling, extracting, rearranging data ...')
                ss = find(ismember(CFG.goodsub, subjIdx)); % which subject fibre extract from allrois{roi}
                roidata = allrois{iroi}(:,:,ss);
                roidata = permute(roidata(:,:,:),[1 3 2]); % trick to do reshape correctly
                roidata = reshape(roidata, CFG.nclust*length(subjIdx), CFG.nfreq); 

                disp('Calculating nlogl ...')
                for jroi = CFG.goodroi
                    X                   = roidata; % nsub*10 x 42
                    obj                 = gm{jroi}; % gaussian model class
                    if(obj.NumComponents == 0)
                        post = 0;
                        nlogl = inf;
                    else
                        % obj - 2nd lvl GM model for ROI, X - given ROI data from all subjects (pooled above) from test set
                        [post, nlogl]       = posterior(obj, X);
                        roifit(iroi, jroi)  = nlogl;
                    end
                end
                
                disp('Calculating ranks ...')
                roifitrank(iroi,CFG.goodroi) = tiedrank(roifit(iroi,CFG.goodroi)); % tiedrank gives rank 1 to the smallest value (so it is perfect for finding minimum negative log likelihood)
            end


        case 'subj_avg'
            disp('SUBJ_AVG mode !!!')
            % CALCULATE POSTERIORS
            disp('CALCULATING POSTERIORS and Negative Log-likelihood. NlogL will be ranked additionaly.')

            % later usage of nanmean will account for missing subjects (so
            % CFG.goodsub can be arbitrary)
            roifit      = nan(CFG.nroi, CFG.nroi, CFG.nsub);
            roifitrank  = nan(CFG.nroi, CFG.nroi, CFG.nsub);

            isub = 1;
            iroi = 1;
            jroi = 1;
            disp('Entering subjects loop ...')
            for isub = 1:length(subjIdx)
                disp(['Analysing subject ...', num2str(subjIdx(isub)) ' ...'])
                
                disp('Entering goodroi loop ...')
                for iroi = CFG.goodroi
                    iroi

                    % pool data
                    % POOLING - rearrange data (all 11 subject clusters in 1st dimension, frequencies in 2nd dimension)
                    % (Pooling is my assumption, because this detail was not provided in the paper.)
                    
                    disp('Pooling, extracting, rearranging data ...')
                    ss = find(ismember(CFG.goodsub, subjIdx)); 
                    roidata = allrois{iroi}(:,:,ss);
                    roidata = permute(roidata(:,:,isub),[1 3 2]); % trick to do reshape correctly
                    roidata = reshape(roidata, CFG.nclust*length(isub), CFG.nfreq); 

                    disp('Calculating nlogl ...')
                    for jroi = CFG.goodroi
                        obj                         = gm{jroi}; % gaussian model class
                        X                           = roidata; 
                        % obj - 2nd lvl GM model for ROI, X - given ROI data from all subjects (pooled above) from test set
                        [post, nlogl]               = posterior(obj, X);

                        roifit(iroi, jroi, isub)    = nlogl;
                    end
                    disp('Calculating ranks ...')
                    roifitrank(iroi, CFG.goodroi, isub) = tiedrank(roifit(iroi, CFG.goodroi, isub)); % tiedrank gives rank 1 to the smallest value (so it is perfect for finding minimum negative log likelihood)
                end

            end
            disp('Averaging across subjects ...')
            % average across subjects
            roifit      = nanmean(roifit,3);
            roifitrank  = nanmean(roifitrank,3);
    end

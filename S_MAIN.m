% S_MAIN
% Started 16.07.2017
% Michal Komorowski, michak@is.umk.pl

% These script available in github repository are meant to be used to reproduce results in paper “Individual Human Brain Areas Can Be Identified from Their Characteristic Spectral Activation Fingerprints”  by A. Keitel & J. Gross 2016 with data publicly available thanks to the authors.

% Assuming that you have MATLAB with Fieldtrip toolbox and read whole script (especially comments) first to know how to configure (CONFIGRURATION SECTION) and prepare directories with data it to run it smoothly.
% If there are any questions, please contact me - michak@is.umk.pl


clear; close all; clc;

% !!! Initialize fieldtrip here !!!

clear; close all; clc;
%% %%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atlas: Note that this is a template from an earlier fieldtrip version, these
% templates change. Depending on your template/spatial resolution, the
% interpolation with the atlas will be different. % -> THIS CAN AFFECT GM MODELS !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%% CONFIGRURATION SECTION %%%%%%%%%%%%%%%%%
global CFG
% -!- Free to edit ---
CFG.HomeFolder          = './';
CFG.rawDir              = [CFG.HomeFolder,'../author_f/'] % <- path to the raw data files (clean MEG, LCMV, templategrid10mm.mat)
CFG.rsltsDir            = [CFG.HomeFolder, 'data/', 'data4a/'] % results Dir !!! Important for calcultation and figures !!!
CFG.FieldtripPath       = '...';  % <- provide path to your fieldtrip directory !!!
CFG.AtlasDir            = [CFG.FieldtripPath 'template/atlas/aal/ROI_MNI_V4.nii']; % templates for atlas (AAL atlas) and 10-mm source model; S_Ranking_Homologue is hard-code-dependent of AAL atlas with 116 ROIs !
CFG.SourceMod           = [CFG.rawDir 'templategrid10mm.mat'];
CFG.goodsub             = 1:22; %1:22; -> all
CFG.goodroi             = [1:94 96:116];   
CFG.nfreq               = 42;  % number of analysed frequencies after fieldtrip takes out the redundant ones (check numel(freq.freq) after FFT)
% GM1
CFG.regu                = 0.012;
CFG.indivname           = 'indivclust';
CFG.nclust              = 10;  % number of individual clusters for 1-st-level models
% GM2
CFG.optname             = sprintf('%soptclust', CFG.rsltsDir);
CFG.clustnum            = 1:15; % list of number of clusters to evaluate (during optimal number of clusters in 2nd lvl GM model)
CFG.niter               = 1000; % !!! DIVISIBLE by 4 !!! -> look at S_OptClustPARA.m and S_SAVEOptClustPARA.m number of iterations for optimal cluster analysis 
CFG.finalname           = sprintf('%sGMmodel_allsubs', CFG.rsltsDir);
CFG.majorN              = 16; % how many subjects has to contribute to 2nd lvl cluster to be accepted (set N=1 for accepting all clusters)
                               % WARNING ! This number should change if you change the number of subjects - check Chi^2 in the paper.
% Classification
CFG.Niter               = 120; % how many iterations in classification
CFG.fixclnum            = 4; % fixed number of cluters in 2nd lvl GM used in classification (4 was the median mentioned in article)
CFG.NlogLRankMode       = 'pool'; % look at description in F_NlogLRank.m function

% -!- do not edit ---
CFG.logfreqs            = logspace(0,log10(120),50); 
CFG.nsub                = length(CFG.goodsub); % opposing to CFG.nroi which is always maximum number of ROI in used atlas, CFG.nsub is the number of accepted subjects; You can get subj index with CFG.goodsub(ss)
ft_defaults             % initialize fieldtrip
addpath(CFG.FieldtripPath);
CFG.atlas               = ft_read_atlas(CFG.AtlasDir);
CFG.atlas               = ft_convert_units(CFG.atlas, 'cm');
CFG.nroi                = numel(CFG.atlas.tissuelabel);        % number of ROIs in atlas
load(CFG.SourceMod); 
cfg = [];
cfg.interpmethod        = 'nearest';
cfg.parameter           = 'tissue';
CFG.sourceatlas         = ft_sourceinterpolate(cfg, CFG.atlas, sourcemodel);
clear sourcemodel cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%% PLOTTING SECTION %%%%%%%%%%%%%%%%%%%%%%%

% !!! Run CONFIGRURATION SECTION first !

S_Fig2          % ROI spectral profiles (spectrum, surf)
S_Fig5          % Distribution of number of clusters (hist, surf, regr)
S_Fig3b         % Classification procedure and results (hist, surf, regr)
S_Fig3b_homo    % Classification procedure and results (hist, surf, regr) - accepting homologues
S_Fig4          % Area networks according to similarity analysis (surf, slice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%% MAIN CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%

%% --- 1st lvl GM ---
% This is LONG step.

for isub = CFG.goodsub
    fprintf('\nAnalysing sub %d .. \n', isub);
    
    %(ff - analysed frequencies by ft_freqanalysis - Fieldtrip could force little bit different freqs than those provided by the user)
    [sourcef, ff] = F_SourceProjNorm(isub); 
                                          
    allsourcef{isub} = sourcef; 
    ff_all{isub} = ff;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [roicluster, pt] = F_GM1(CFG.goodroi, sourcef); 
    save(sprintf('%s%s_nclust%d_sub%d', CFG.rsltsDir, CFG.indivname, CFG.nclust, isub), 'roicluster','pt'); % !!! name 'nclust' so files for F_Rearrange('main')
    
end
save([CFG.rsltsDir, 'allsourcef'], 'allsourcef', '-v7.3') % about 8 GB !!!
save([CFG.rsltsDir, 'ff_all'], 'ff_all') % for later plotting ff for individual
%---------------------------------

%% --- 2nd lvl GM ---
% You can start here with only CFG in workspace
[allrois, subindex] = F_Rearrange('main'); % respawn point!
size(allrois), size(subindex)

% This is VERY VERY LONG step (do it using 4-6 MATLABS).
S_OptClustPARA % do it 4 times, then skip it (YOU have to change iteration indices inside to run it !)
S_SAVEOptClustPARA % do it only once, then skip it

% This is not very long step.
[gm, nsubcl, allidx] = F_GM2(allrois, 'opt', 0, CFG.goodsub, subindex)


% Here we accept cluster present in majority of participants.
gm2 = F_AcceptClusters(gm, nsubcl, CFG.majorN);
gm = gm2;

% CHECK: if any roi has zero gm model with zero components (that could
% happen and it is ok ... )
gm_ok = zeros(1, length(gm));
for iroi = CFG.goodroi
    if(gm{iroi}.NumComponents > 0 )
        gm_ok(iroi) = 1;
    end
end

save(CFG.finalname, 'gm','nsubcl', 'allidx', 'subindex', 'gm_ok')
%---------------------------------

%% --- CLASSIFICATION ---
%clear; close all; clc;

% load CFG again ! (manually)

% UPDATE 20.07.2017 - According to paper they do not calculate 1GM again -
% but this is not a problem (1GM models should not differ much if you
% calculate them again). You can thus substitute data 'nclust' to
% 'testdata'.
% economic version - compute 1GM once for all (full version would be in loop below, calculating 1GM for each partition separately)

for isub = CFG.goodsub
    fprintf('\nAnalysing sub %d .. \n', isub);
    
    
    [sourcef ~] = F_SourceProjNorm(isub);
    
    allsourcef{isub} = sourcef; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [roicluster, pt] = F_GM1(CFG.goodroi, sourcef); 
    save(sprintf('%s%s_testdata%d_sub%d', CFG.rsltsDir, CFG.indivname, CFG.nclust, isub), 'roicluster','pt'); % !!! notice name 'testdata', so it loads files for F_Rearrange('class')
end
save([CFG.rsltsDir, 'allsourcef_class'], 'allsourcef', '-v7.3') % about 8 GB !!!


% respawn point ! You can start here ! RELOAD CFG first !!!
[allrois, subindex] = F_Rearrange('main'); % changing to 'main' will load data computed before classification section 
% This is LONG step (30 mins).
for jj = 1:CFG.Niter
    jj
    permIdx = randperm(CFG.nsub)
    tr_idx = CFG.goodsub(permIdx(1:CFG.nsub/2))
    te_idx = CFG.goodsub(permIdx((CFG.nsub/2 + 1):CFG.nsub))
    
    [gm, ~, ~] = F_GM2(allrois, 'fix', CFG.fixclnum, tr_idx, subindex)
    [roifit, roifitrank] = F_NlogLRank(allrois, gm, te_idx, CFG.NlogLRankMode); 
    iterdata(jj).roifit = roifit;
    iterdata(jj).roifitrank = roifitrank;
end
S_Ranking % get: 1)ranks: ROI match (successful recognition) rank in all iterations  ; compute: 2) ROI match mean rank across all areas and iterations 
S_Ranking_Homologue % as above, but homologue match is considered as successful recognition

save([CFG.rsltsDir, 'classresults'], 'iterdata', 'gm', 'te_idx', 'tr_idx', 'ranks', 't20meanrank', 'ranks_homo', 't20meanrank_homo')
%---------------------------------


%% --- DENDROGRAM ---

% respawn point !
% reload CFG
% load 2nd GM models
load(CFG.finalname)

[allrois, subindex] = F_Rearrange('main'); % it should be set to 'main' if there was 'main' before. Default is 'class'.
% perform NlogL on all subject data
[roifit, ~] = F_NlogLRank(allrois, gm, CFG.goodsub, 'pool'); % should stay with 'pool' according to paper (methods section: "Similarity across Areas")
S_ROIsimil
save([CFG.rsltsDir,'ROIsimilresults'], 'roifit', 'Z', 'T', 'rows')

%---------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%% END ROUTINES %%%%%%%%%%%%%%%%%%%%%%
% save config file if you wish:
choice = questdlg('Do you want to save CFG variable?', '_', 'Yes', 'No', 'No');
switch choice
    case 'Yes'
        save([CFG.rsltsDir,'CFG'], 'CFG')
        disp('Config variable saved.')
    case 'No'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

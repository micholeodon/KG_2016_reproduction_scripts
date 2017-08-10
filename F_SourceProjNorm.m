function [sourcef, ff] = F_SourceProjNorm(isub)
global CFG
    % load clean, segmented data
    %load(sprintf('%sdata_clean_%d',HomeFolder,isub));     % [data]
    load(sprintf('%sdata_clean_%d', CFG.rawDir, isub));     % [data]

    % FFT
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'fourier';
    cfg.taper       = 'dpss';
    cfg.tapsmofrq   = 2; % the amount of spectral smoothing through multi-tapering (+/- df)
    cfg.pad         = 2; % length in seconds to which the data can be padded out.
    cfg.keeptrials  = 'yes';
    cfg.foi         = CFG.logfreqs;
    freq = ft_freqanalysis(cfg, data);

    ff = freq.freq;                                 % save frequency resolution later
    % load previously computed LCMV filter
    load(sprintf('%sLCMV_%d', CFG.rawDir, isub),'filt');     % [filt]
    
    ninsvox = size(filt,1)                        % number of voxels inside brain

    % facilitate data handling by numerically increasing values
    % (can be different for each MEG system; look at data)
    freq.fourierspctrm = freq.fourierspctrm*1e14;

    
    sourcef = zeros(length(freq.cumsumcnt), ninsvox, CFG.nfreq);    % create source matrix (trials x voxels x freq)

    disp('Source projection ...')
    % SOURCE PROJECTION; apply weights to fourier spectra and compute power
    m = 1;
    for k = 1:length(freq.cumtapcnt),
        for k2 = 1:freq.cumtapcnt(1),
            
            
            sourcef(k,:,:) = squeeze(sourcef(k,:,:)) + abs(filt*squeeze(freq.fourierspctrm(m,:,:))).^2;
            m = m+1;
        end
    end

    disp('Normalize activity ...')
    % NORMALISE ACTIVITY (we use ratio normalisation here)
    avgact = squeeze(mean(mean(sourcef,1),2));      % average activity of whole brain
    tmp = repmat(avgact',[ninsvox 1]);
   
    for k = 1:length(freq.cumsumcnt),        
        sourcef(k,:,:) = (squeeze(sourcef(k,:,:)))./tmp;
    end

disp('Source projection DONE !')

end

%

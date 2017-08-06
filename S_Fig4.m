% S_Fig4      % Area networks according to similarity analysis (surf, slice)
addpath ../misc % add some tools
mkdir([CFG.rsltsDir, 'fig/'])

% prepare colormap
disp('Preparing colormap ...')
hexcol = {'#000000', '#00FF00', '#0000FF', '#FF0000', '#01FFFE', '#FFA6FE', '#FFDB66', '#006401', '#010067', '#95003A', '#007DB5', '#FF00F6', '#FFEEE8', '#774D00', '#90FB92', '#0076FF', '#D5FF00', '#FF937E', '#6A826C', '#FF029D', '#FE8900', '#7A4782', '#7E2DD2', '#85A900', '#FF0056', '#A42400', '#00AE7E', '#683D3B', '#BDC6FF', '#263400', '#BDD393', '#00B917', '#9E008E', '#001544', '#C28C9F', '#FF74A3', '#01D0FF', '#004754', '#E56FFE', '#788231', '#0E4CA1', '#91D0CB', '#BE9970', '#968AE8', '#BB8800', '#43002C', '#DEFF74', '#00FFC6', '#FFE502', '#620E00', '#008F9C', '#98FF52', '#7544B1', '#B500FF', '#00FF78', '#FF6E41', '#005F39', '#6B6882', '#5FAD4E', '#A75740', '#A5FFD2', '#FFB167', '#009BFF', '#E85EBE'};
rgbcol = hex2rgb(hexcol);
rgbcol(1,:) = []; % discard black colour
mycolormap = rgbcol;

disp('Loading results ...')
load([CFG.rsltsDir, 'ROIsimilresults'])

disp('Reading MRI ...')
mrifile = [CFG.FieldtripPath 'template/anatomy/single_subj_T1.nii']
mri = ft_read_mri(mrifile)
mri.coordsys = 'mni'; % to prevent manual fixing of coordsys
%ft_sourceplot([],mri) % uncomment if you want to see raw mri in ortho view

disp('Interpolate to atlas ...')
cfg = [];
cfg.parameter = 'tissue';
interp = ft_sourceinterpolate(cfg, CFG.atlas, mri)

% Because T vector is size of length(CFG.goodroi) because not all rois are good, so indexing is disturbed.
% We need to restore correct indexing inserting missing rows in T with NaNs.
disp('Correcting indexing ...')
Tcorr = nan(1,CFG.nroi)
Tcorr(CFG.goodroi) = T


% colour voxel by their t20meanROIrank
disp('Coloring ROIs ...')
interp.rank = interp.tissue
for iroi = CFG.goodroi
    roivoxs = find(interp.tissue==iroi);
    interp.rank(roivoxs) = Tcorr(iroi);       
end

disp('Creating surface view ...')
cfg                     = [];
cfg.method              = 'surface';
cfg.projmethod          = 'project';
cfg.camlight            = 'yes';
cfg.locationcoordinates = 'voxel';
cfg.cmap                = mycolormap(1:20,:); % color according to number of ROI
cfg.cmap                = [[0,0,0]; cfg.cmap] % color map
cfg.funcolormap         = cfg.cmap;
cfg.funparameter        = 'rank';
cfg.funcolorlim         = [0 20];
cfg.atlas               = CFG.atlas;
ft_sourceplot(cfg, interp) % weird black spots - to be solved soon
title('Area networks')
saveas(gcf,[CFG.rsltsDir, 'fig/an_surf_top', '.png'])
view([90 0])
saveas(gcf,[CFG.rsltsDir, 'fig/an_surf_right', '.png'])
view([-90 0])
camlight
saveas(gcf,[CFG.rsltsDir, 'fig/an_surf_left', '.png'])


% 1) Slice
disp('Creating slice view ...')
cfg                     = [];
cfg.method              = 'slice';
cfg.nslices             = 20;
cfg.slicerange          = [1 81];
cfg.locationcoordinates = 'voxel';
cfg.location            = [41,52,49]; % some custom typical location (in voxel, not cm) !
cfg.cmap                = mycolormap(1:20,:); % color according to number of ROI
cfg.cmap                = [[0,0,0]; cfg.cmap] % color map
cfg.funcolormap         = cfg.cmap;
cfg.funparameter        = 'rank';
cfg.funcolorlim         = [0 20];
cfg.atlas               = CFG.atlas;
ft_sourceplot(cfg, interp)
title('Area networks')
saveas(gcf,[CFG.rsltsDir, 'fig/an_slice', '.png'])



% -----------------------------------------------
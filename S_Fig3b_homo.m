% S_Fig3b but when considering homologue areas as correctly recognized
% area.

mkdir([CFG.rsltsDir, 'fig/'])

%% %%% MAIN CODE %%%

%% --- Regression plot ---
% distance from brain "center" (from Supplementary Materials)
r_ROI = [ 6.280083847 6.689159586 5.4561476634 5.5521367675 5.7445626465 5.5277079839 5.6910814706 5.9760459336 5.5897135847 5.9710184411 5.2586178055 5.4637048791 5.5445014743 6.0004616628 4.8455305867 5.3346352578 5.0159744816 5.3397000566 5.8941767402 6.1536450124 1.8708286934 2.7726341266 5.5641120191 5.7542255005 5.123475383 5.1613951602 4.1570722871 4.1533119315 3.6929331299 3.8703477669 3.4456604849 4.2941821107 4.6720578978 3.8649421442 5.2230259429 4.9117206761 3.8262252939 3.9201207811 3.4776069358 3.7027017163 3.2015621187 3.6055512755 8.1905433275 7.7194739314 8.3784453873 8.402239769 6.6964047619 7.1340774967 9.3282650544 9.2048900048 8.9502522567 9.1292146734 8.7401730532 9.1782351245 5.3233239407 5.6286268701 6.9270794037 7.1505520666 8.7611536055 8.89756521 7.7365973316 7.874007874 7.4393380082 7.4458041876 8.727705283 8.4440789395 7.3196862038 7.1111996877 7.4402973523 7.4760032772 1.8126539343 1.9309052441 2.1852940773 2.7888667551 1.5 2 2.536158269 2.487003254 5.2704627669 5.0249378106 5.6400354609 6.3403686287 4.6583258795 5.3062463192 6.6337557142 6.9575043173 5.3471335063 5.7183913822 6.3393611666 6.78603683 7.9916956899 8.377346286 8.8789388755 8.7559979443 NaN 4.582575695 4.9764879281 4.9843505093 6.7750148669 7.0953329732 7.9882726543 9.5 7.8164712129 8.1965775966 7.0710678119 7.1721684308 5.3851648071 5.8309518948 4.472135955 4.1231056256 5.4589376256 6.8718427094 7.9056941504 7.3824115301 6.8007352544 5.8309518948]';


load([CFG.rsltsDir, 'classresults'])

t20meanROIrank_homo = trimmean(ranks_homo,20,2);
% put data in the form of regression
Y = t20meanROIrank_homo(CFG.goodroi);
X = [ ones(length(r_ROI),1) r_ROI ];
X = X(CFG.goodroi,:);


% (p - number of coefficients, here 2; n - number of data points)
% bint - p-by-2 matrix bint of 95% confidence intervals for the coefficient estimates
% r    -  n-by-1 vector r of residuals.
% rint - n-by-2 matrix rint of intervals that can be used to diagnose outliers.
% stats - = [R2 F pF evar] R^2 statistic, the F statistic and its p value, and an estimate of the error variance.
% !!! Assumption - stats are calculated for model with intercept term !
[b, bint, r, rint, stats] = regress(Y,X)
R2 = stats(1)
F = stats(2)
pF = stats(3) 
evar = stats(4)

Yfit = X*b;

fig = figure
scatter(X(:,2), Y)
hold on 
plot(X(:,2), Yfit)
hold off
xlabel('Radius r (in cm)')
ylabel('Mean rank / area')
set(gca,'FontSize',18)
msg = ['y = (', num2str(b(1)), ') + (' , num2str(b(2)), ') * x']


title(msg)
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_r_corr', '.png'])

% -----------------------------------------------

%% --- Histogram ---

his = histogram(round(t20meanROIrank_homo),'BinMethod', 'integers')
his = his.Values;

N = length(his);

colors = parula(100);
close
figure;
for i=1:N
    bb = bar(i, his(i), 1);
    bb.FaceColor = colors(i*floor(100/N),:);
    hold on;
end
set(gca,'Xtick',1:N,'XTickLabel', 1:N)
xlabel('Mean rank in classification')
ylabel('Number of areas')
set(gca,'FontSize',14)
set(gcf, 'Position', [0 0 1200 800])
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_hist', '.png'])
% -----------------------------------------------

%% --- Brain surface ---

mrifile = [CFG.FieldtripPath 'template/anatomy/single_subj_T1.nii']
mri = ft_read_mri(mrifile)
mri.coordsys = 'mni'; % to prevent manual fixing of coordsys
%ft_sourceplot([],mri)

cfg = [];
cfg.parameter = 'tissue';
interp = ft_sourceinterpolate(cfg, CFG.atlas, mri)

% colour voxel by their t20meanROIrank
interp.rank = interp.tissue
for iroi = 1:CFG.nroi
    roivoxs = find(interp.tissue==iroi);
    interp.rank(roivoxs) = t20meanROIrank_homo(iroi);       
    t20meanROIrank_homo(iroi)
end


cfg                     = [];
cfg.method              = 'surface';
cfg.projmethod          = 'project';
cfg.camlight            = 'yes';
cfg.locationcoordinates = 'voxel';
cfg.cmap                = parula; % color according to number of ROI
cfg.cmap                = [[0,0,0]; cfg.cmap] % color map
cfg.funcolormap         = cfg.cmap;
cfg.funparameter        = 'rank';
cfg.funcolorlim         = [0 4];
cfg.atlas               = CFG.atlas;

cfg.surffile = [CFG.FieldtripPath 'template/anatomy/surface_pial_both.mat']; % weird black spots - to be solved soon
ft_sourceplot(cfg, interp)
camlight('right')
title('Mean rank')
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_both', '.png'])

cfg.surffile = [CFG.FieldtripPath 'template/anatomy/surface_pial_left.mat']; % weird black spots - to be solved soon
ft_sourceplot(cfg, interp)
view([90 0])
camlight('right')
title('Mean rank')
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_left_in', '.png'])

ft_sourceplot(cfg, interp)
view([-90 0])
camlight('right')
title('Mean rank')
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_left_out', '.png'])

cfg.surffile = [CFG.FieldtripPath 'template/anatomy/surface_pial_right.mat']; % weird black spots - to be solved soon
ft_sourceplot(cfg, interp)
view([90 0])
camlight('right')
title('Mean rank')
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_right_out', '.png'])

ft_sourceplot(cfg, interp)
view([-90 0])
camlight('right')
title('Mean rank')
saveas(gcf,[CFG.rsltsDir, 'fig/mrank_homo_right_in', '.png'])

% -----------------------------------------------
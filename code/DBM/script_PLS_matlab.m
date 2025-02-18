% %script_PLS.m
% %Setup a parallel pool to make matlab faster in the CIC cluster
% number_of_cores=12;
% d=tempname(); %get a temporary location;
% mkdir(d);
% cluster = parallel.cluster.Local('JobStorageLocation',d,'NumWorkers',number_of_cores);
% parpool(cluster,number_of_cores);

% This script uses the Partial Least Squares from Baycrest, documentation> https://www.rotman-baycrest.on.ca/index.php?section=100

%% Setting up
cd '/data/chamal/projects/jalilr/pls'
addpath /data/chamal/projects/jalilr/pls % point to path to directory where you will run the PLS
addpath /data/chamal/projects/jalilr/pls/plscmd %point to path with PLS tools
inpath = '/data/chamal/projects/jalilr/pls/mice_data/relative_jacobians/nii/t3_pls/' % point to path to directory where you have the niftis that you will run the PLS on
data_dir = '/data/chamal/projects/jalilr/pls/outputs';

%% Import data
data = readtable(fullfile('/data/chamal/projects/jalilr/pls/DLC_tracked_4pls_path.csv'));

%for behavs you can change the variable's names
behavs = {'NOR_PI','NOR_PI_entries','NOR_Interaction_index','NOR_nose.distance.moving','NOR_nose.speed.moving','NOR_nose.time.moving','NOR_nose.time.stationary','NOR_nose.ob.distance.moving','NOR_nose.ob.speed.moving','NOR_nose.ob.time.moving','NOR_nose.ob.total.time','NOR_nose.ob.time.stationary','NOR_nose.ob.transitions','NOR_nose.nob.distance.moving','NOR_nose.nob.speed.moving','NOR_nose.nob.time.moving','NOR_nose.nob.total.time','NOR_nose.nob.time.stationary','NOR_nose.nob.transitions','EPM_PI','EPM_PI_entries','EPM_Anxiety_index','EPM_nose.dip','EPM_bodycentre.distance.moving','EPM_bodycentre.speed.moving','EPM_bodycentre.time.moving','EPM_bodycentre.time.stationary','EPM_bodycentre.center.distance.moving','EPM_bodycentre.center.speed.moving','EPM_bodycentre.center.time.moving','EPM_bodycentre.center.total.time','EPM_bodycentre.center.time.stationary','EPM_bodycentre.center.transitions'}
nbehav = numel(behavs); % gets the number of variables that you are interested in (those specified above)
path_nii = unique(data.relative_jacobian); %path to nii files


%% MANUALLY IMPORT MATLAB ANALYSIS/ MATLAB DATA
% load('/data/chamal/projects/jalilr/pls/outputs/PLS.mat')

%% Creating variables
allsubjects = unique(data.RID); % get subject ids
naidx = strcmp('NA',data.relative_jacobian); % find subjects that don't have jac. file
nasubj = data.RID(naidx); % count the number of subj that don't have jac file
subjects = setxor(nasubj,allsubjects); % make array with only the subjects that have a jac file
groups = unique(data.Group); %make array with the group names

%% Load the DBM brain mask
mask_nii = niftiinfo(fullfile(inpath,'mask.nii.gz'));  % must in same directory as nifit-jacobians
mask = niftiread(mask_nii);
mask_idx = find(mask);
mri_path = fullfile(data.relative_jacobian);

%% Setting up MRI variables
ngroup = numel(groups); % get the number of groups
nsubj = numel(subjects); % get the number of subjects
nvoxel = length(mask_idx); % get the number of DBM measures (ie the number of voxels)

%% Creating matrix with the imaging data
mri_data = zeros(nvoxel, nsubj); %fill with zeros a matrix of subjects times mask's number of voxels

%%for each subject gather nifti path and gather the info of each voxel within the mask
for i=1:nsubj
 nii= niftiinfo(mri_path{i,:});
 niidata = niftiread(nii);
 nifti_mask = niidata(mask_idx);
 mri_data(:,i)=nifti_mask;
end

% defining neuroimaging data normalized
mri_data_ok = mri_data.';


%% setting up behaviour variables

%selecting the behaviours from column 2 to 7, and from column 9 to 47
data_beh = table2array(data(:,[3:35]));

%%Normalizing behavioural variables
mub1=mean(data_beh);
sigmab1 = std(data_beh);
zb1 = bsxfun(@minus,data_beh,mub1);
zscored_behavdata = bsxfun(@rdivide,zb1,sigmab1);

% defining behavioral data normalized
behav = zscored_behavdata(1:nsubj,1:nbehav);

%% run PLS

% matrix with behavioral data
behav_data{1} = behav;
data_behav{1} = [behav_data{1}];

%matrix with neuroimaging data
mri_data_1{1} = mri_data_ok;
datamat{1} = [mri_data_1{1}];


% Setting PLS parameters

option.method = 3;  % 1 = mean-centering (i.e. group/condition comparison)
                    % 3 = behavioural PLS (comparing 2 sets of variables)

option.num_perm = 1000;     %number of permutations
% option.num_split = 100;
option.num_boot = 1000;     %number of bootstraps

%save behav_data in the stacked matrix
option.stacked_behavdata = [behav_data{1}];

result = pls_analysis(datamat,nsubj,1,option); % (dataframe, #subj,#condition, option)


%save matlab analysis, change to correct path
save('/data/chamal/projects/jalilr/pls/outputs/PLS.mat');


%% check p-values and % covariance

%check p-values
result.perm_result.sprob
% ans =

    0.1668
    0.3836
    0.2887
    0.0160 *
    0.7033
    0.9241
    0.9720
    0.6893
    0.3387
    0.9590
    0.9950
    0.8561
    0.9890
    0.8611
    0.9980
    0.9990
    0.9750
    0.9940
    0.9920
    0.9990
    0.9371


% percent covariance
pct_cov = (result.s .^2) / sum(result.s .^2);
pct_cov

    0.3117
    0.2307
    0.1452
    0.0897 *
    0.0581
    0.0391
    0.0299
    0.0257
    0.0203
    0.0152
    0.0112
    0.0089
    0.0058
    0.0040
    0.0025
    0.0012
    0.0005
    0.0004
    0.0001
    0.0000
    0.0000

% result.perm_splithalf.ucorr_prob
% result.perm_splithalf.vcorr_prob

%% Figures




% visualize the %covariance for each LV

figure;
% Set default font
set(gca, 'FontName', 'Helvetica')
% Left y-axis
yyaxis left
plot(1:20, pct_cov, 'k.', 'MarkerSize', 20, 'LineWidth', 2) % numbers are referring to LVs
xlabel('Latent Variables', 'FontSize', 15, 'FontName', 'Helvetica')
ylabel('Percent Variance Explained', 'FontSize', 15, 'FontName', 'Helvetica')
set(gca, 'YColor', 'k') % Set left y-axis color to black
% Right y-axis
yyaxis right
color_right =  [214/255, 39/255, 40/255];
plot(1:20, result.perm_result.sprob, '.', 'Color', color_right, 'MarkerSize', 30, 'LineWidth', 2)
ylabel('p-value', 'FontSize', 15, 'FontName', 'Helvetica')
ylim([0 0.05])
xlim([0 20])
set(gca, 'YColor', color_right) % Set right y-axis color to black
% Increase overall font size
set(gca, 'FontSize', 15, 'FontName', 'Helvetica')
exportgraphics(gcf,'cov_vs_pval.png','Resolution',300)



% plot split-half resampling probability for each LVs

% figure;
% yyaxis left
% plot(1:45,result.perm_splithalf.ucorr_prob, 'b.', 'MarkerSize',15,'LineWidth',2) %numbers are referring to LVs
% xlabel('Latent Variables','FontSize',15)
% ylabel('P U correlation','FontSize',15)
% ylim([0 0.08])
% yyaxis right
% plot(1:24,result.perm_splithalf.vcorr_prob, '.', 'MarkerSize',15,'LineWidth',2)
% ylabel('P V correlation','FontSize',15)
% ylim([0 0.08])
% xlim([0 45])
% exportgraphics(gcf,'split_prob.png','Resolution',300)





% LATENT VARIABLES FIGURES********************
% NOTE
% You can choose among different measures of association between brain data and behavioral data, Pearson correlation is the default. X-axis reflects the correlations between behavior data and brain scores. If bootstrap test is included, correlation bar will come with an error bar specified by upper and lower error range for the correlation.


names = {'sex_m','NOR_PI','NOR_PI_entries','NOR_Interaction_index','NOR_nose.ob.total.time','NOR_nose.ob.time.stationary','NOR_nose.ob.percentage.moving','NOR_nose.ob.transitions','NOR_nose.nob.time.moving','NOR_nose.nob.total.time','NOR_nose.nob.time.stationary','EPM_bodycentre.center.distance.moving','EPM_bodycentre.center.time.moving','EPM_bodycentre.center.percentage.moving','EPM_bodycentre.closed.raw.speed','EPM_bodycentre.closed.percentage.moving','EPM_bodycentre.closed.top.raw.speed','EPM_bodycentre.closed.top.percentage.moving','EPM_bodycentre.closed.bottom.raw.speed','EPM_bodycentre.closed.bottom.percentage.moving','EPM_bodycentre.open.right.speed.moving'}
names = behavs;
color_demo = [31/255, 119/255, 180/255];
color_behavtest1 =  [255/255, 219/255, 88/255];
color_behavtest2 = [255/255, 127/255, 14/255];
color_gray = [137 137 137] / 255;
width =0.9; % Adjust this value as needed

figure;
demo = result.lvcorrs(1,1);
behavtest1 = result.lvcorrs(2:5,1);
behavtest2 = result.lvcorrs(6:50,1);
upper = result.boot_result.ulcorr(:,1) - result.lvcorrs(:,1);
lower = result.lvcorrs(:,1) - result.boot_result.llcorr(:,1);
barh(1,demo, 'FaceColor', color_demo, 'BarWidth', width); hold on
barh(2:5,behavtest1,'FaceColor', color_behavtest1, 'BarWidth', width); hold on
barh(6:50,behavtest2,'FaceColor', color_behavtest2, 'BarWidth', width); hold on
set(gca,'TickLabelInterpreter','none')
set(gca,'Ytick',1:numel(names),'YTickLabel',names,'FontSize',8,'FontName', 'Helvetica');
er = errorbar(result.lvcorrs(:,1),1:length(result.lvcorrs(:,1)),lower,upper,'horizontal');
er.LineStyle = 'none';
er.Color = 'k';
set(gca, 'TickLabelInterpreter', 'none', 'YTick', 1:numel(names),
    'YTickLabel', names, 'FontSize', 14, 'FontName', 'Helvetica',
    'Box', 'off', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5);
title ('LV1 with error bars from bootstrap resampling','FontSize', 13, 'FontName', 'Helvetica')
% Customize grid and axes
grid on;
ax = gca;
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0.1 0.1 0.1]; % Dark gray grid lines
ax.GridAlpha = 0.2; % Grid line transparency
% Remove background color (like `theme_bw()`)
set(gca, 'Color', 'w');
xlim([-1 1])
hold off
exportgraphics(gcf,'LV1_pearsoncorr.png','Resolution',300)


figure;
demo = result.lvcorrs(1,2);
behavtest1 = result.lvcorrs(2:5,2);
behavtest2 = result.lvcorrs(6:50,2);
upper = result.boot_result.ulcorr(:,2) - result.lvcorrs(:,2);
lower = result.lvcorrs(:,2) - result.boot_result.llcorr(:,2);
barh(1,demo, 'FaceColor', color_demo, 'BarWidth', width); hold on
barh(2:5,behavtest1,'FaceColor', color_behavtest1, 'BarWidth', width); hold on
barh(6:50,behavtest2,'FaceColor', color_behavtest2, 'BarWidth', width); hold on
set(gca,'TickLabelInterpreter','none')
set(gca,'Ytick',1:numel(names),'YTickLabel',names,'FontSize',8,'FontName', 'Helvetica');
er = errorbar(result.lvcorrs(:,2),1:length(result.lvcorrs(:,2)),lower,upper,'horizontal');
er.LineStyle = 'none';
er.Color = 'k';
set(gca, 'TickLabelInterpreter', 'none', 'YTick', 1:numel(names),
    'YTickLabel', names, 'FontSize', 14, 'FontName', 'Helvetica',
    'Box', 'off', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5);
title ('LV2 with error bars from bootstrap resampling','FontSize', 13, 'FontName', 'Helvetica')
% Customize grid and axes
grid on;
ax = gca;
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0.1 0.1 0.1]; % Dark gray grid lines
ax.GridAlpha = 0.2; % Grid line transparency
% Remove background color (like `theme_bw()`)
set(gca, 'Color', 'w');
xlim([-1 1])
hold off
exportgraphics(gcf,'LV2_pearsoncorr.png','Resolution',300)


% figure;
% demo = result.lvcorrs(1,3);
% behavtest1 = result.lvcorrs(2:5,3);
% behavtest2 = result.lvcorrs(6:50,3);
% upper = result.boot_result.ulcorr(:,3) - result.lvcorrs(:,3);
% lower = result.lvcorrs(:,3) - result.boot_result.llcorr(:,3);
% barh(1,demo, 'FaceColor', color_demo, 'BarWidth', width); hold on
% barh(2:5,behavtest1,'FaceColor', color_behavtest1, 'BarWidth', width); hold on
% barh(6:50,behavtest2,'FaceColor', color_behavtest2, 'BarWidth', width); hold on
% set(gca,'TickLabelInterpreter','none')
% set(gca,'Ytick',1:numel(names),'YTickLabel',names,'FontSize',8,'FontName', 'Helvetica');
% er = errorbar(result.lvcorrs(:,3),1:length(result.lvcorrs(:,3)),lower,upper,'horizontal');
% er.LineStyle = 'none';
% er.Color = 'k';
% set(gca, 'TickLabelInterpreter', 'none', 'YTick', 1:numel(names),
%     'YTickLabel', names, 'FontSize', 14, 'FontName', 'Helvetica',
%     'Box', 'off', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5);
% title ('LV3 with error bars from bootstrap resampling','FontSize', 13, 'FontName', 'Helvetica')
% % Customize grid and axes
% grid on;
% ax = gca;
% ax.GridLineStyle = '-'; % Solid grid lines
% ax.GridColor = [0.1 0.1 0.1]; % Dark gray grid lines
% ax.GridAlpha = 0.2; % Grid line transparency
% % Remove background color (like `theme_bw()`)
% set(gca, 'Color', 'w');
% xlim([-1 1])
% hold off
% exportgraphics(gcf,'LV3_pearsoncorr.png','Resolution',300)


figure;
demo = result.lvcorrs(1,4);
behavtest1 = result.lvcorrs(2:5,4);
behavtest2 = result.lvcorrs(6:50,4);
upper = result.boot_result.ulcorr(:,4) - result.lvcorrs(:,4);
lower = result.lvcorrs(:,4) - result.boot_result.llcorr(:,4);
barh(1,demo, 'FaceColor', color_demo, 'BarWidth', width); hold on
barh(2:5,behavtest1,'FaceColor', color_behavtest1, 'BarWidth', width); hold on
barh(6:50,behavtest2,'FaceColor', color_behavtest2, 'BarWidth', width); hold on
set(gca,'TickLabelInterpreter','none')
set(gca,'Ytick',1:numel(names),'YTickLabel',names,'FontSize',8,'FontName', 'Helvetica');
er = errorbar(result.lvcorrs(:,4),1:length(result.lvcorrs(:,4)),lower,upper,'horizontal');
er.LineStyle = 'none';
er.Color = 'k';
set(gca, 'TickLabelInterpreter', 'none', 'YTick', 1:numel(names),
    'YTickLabel', names, 'FontSize', 14, 'FontName', 'Helvetica',
    'Box', 'off', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5);
title ('LV4 with error bars from bootstrap resampling','FontSize', 13, 'FontName', 'Helvetica')
% Customize grid and axes
grid on;
ax = gca;
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0.1 0.1 0.1]; % Dark gray grid lines
ax.GridAlpha = 0.2; % Grid line transparency
% Remove background color (like `theme_bw()`)
set(gca, 'Color', 'w');
xlim([-1 1])
hold off
exportgraphics(gcf,'LV4_pearsoncorr.png','Resolution',300)



%% plot bootstrapping histograms
% figure;
% histogram(result.boot_result.compare_u(:,1));
% exportgraphics(gcf,'boot_histoLV1.png','Resolution',300)
figure;
histogram(result.boot_result.compare_u(:,2));
exportgraphics(gcf,'boot_histoLV2.png','Resolution',300)
% figure;
% histogram(result.boot_result.compare_u(:,3));
% exportgraphics(gcf,'boot_histoLV3.png','Resolution',300)
figure;
histogram(result.boot_result.compare_u(:,4));
exportgraphics(gcf,'boot_histoLV4.png','Resolution',300)


%% get bootstrap ratio for each LV
% For brain variables - a bootstrap ratio is calculated as the ratio of the singular vector weight of a
% brain var (ie brain score) to its standard error across bootstraps and we usually look at the 95%
% confidence interval, or like 2 standard deviations away (i.e. the tail ends)

%for LV1
bsr = result.boot_result.compare_u(:,1);
bsr_volume = mask;
bsr_volume(mask_idx) = result.boot_result.compare_u(:,1);
niftiwrite(bsr_volume, 'LV1_scores.nii', mask_nii);

% for LV2
bsr = result.boot_result.compare_u(:,2);
bsr_volume = mask;
bsr_volume(mask_idx) = result.boot_result.compare_u(:,2);
niftiwrite(bsr_volume, 'LV2_scores.nii', mask_nii);

% for LV3
bsr = result.boot_result.compare_u(:,3);
bsr_volume = mask;
bsr_volume(mask_idx) = result.boot_result.compare_u(:,3);
niftiwrite(bsr_volume, 'LV3_scores.nii', mask_nii);

% for LV4
bsr = result.boot_result.compare_u(:,4);
bsr_volume = mask;
bsr_volume(mask_idx) = result.boot_result.compare_u(:,4);
niftiwrite(bsr_volume, 'LV4_micedata.nii', mask_nii);

%% get brain and behavior scores for each LV
brainsc_LV1 = datamat{1} * result.u(:,1);
behavsc_LV1 = data_behav{1} * result.v(:,1);

brainsc_LV2 = datamat{1} * result.u(:,2);
behavsc_LV2 = data_behav{1} * result.v(:,2);

brainsc_LV3 = datamat{1} * result.u(:,3);
behavsc_LV3 = data_behav{1} * result.v(:,3);

brainsc_LV4 = datamat{1} * result.u(:,4);
behavsc_LV4 = data_behav{1} * result.v(:,4);

%% Write brain/behav scores to csv

%LV1
writematrix(brainsc_LV1, "brainsc_micedataLV1.csv");
writematrix(behavsc_LV1, "behavsc_micedataLV1.csv");

%LV2
writematrix(brainsc_LV2, "brainsc_micedataLV2.csv");
writematrix(behavsc_LV2, "behavsc_micedataLV2.csv");

% %LV3
writematrix(brainsc_LV3, "brainsc_micedataLV3.csv");
writematrix(behavsc_LV3, "behavsc_micedataLV3.csv");

%LV4
writematrix(brainsc_LV4, "brainsc_micedataLV4.csv");
writematrix(behavsc_LV4, "behavsc_micedataLV4.csv");
%% Extract temporal components from spatial map time series
%	This script extracts the temporally recurrent network states from
% spatial or 


%%	SET UP FILE SYSTEM

% Clear workspace
clear; close all; clc

% Shuffle random seed.  Necessary to avoid repeating random seeds across parallel computations.
rng("default");
rng("shuffle");

% Find general path (enclosing folder of current directory)
pth{1} = string(strsplit(pwd, filesep));
pth{3,1} = fullfile(pth{1}{:});
pth{2,1} = fullfile(pth{1}{1:end-1});
pth{1,1} = fullfile(pth{1}{1:end-2});

% Set data-specific subdirectories
pth{4,1} = fullfile(pth{2}, "Data");
pth{5,1} = fullfile(pth{3}, "Results");

% List relevant paths
fpth{1,1} = "Functions";
fpth{2,1} = fullfile("MATLAB","spm12");
fpth{3,1} = fullfile("MATLAB","dmd-neuro");
fpth{4,1} = fullfile("MATLAB","gift","GroupICAT","icatb");
fpth{5,1} = fullfile("MATLAB","permutationTest");
fpth{6,1} = fullfile("MATLAB","BCT");
fpth{7,1} = fullfile("MATLAB","DataViz","daviolinplot");

% Add relevant paths
addpath(genpath(fullfile(pth{2}, fpth{1})));
addpath(fullfile(pth{1}, fpth{2}));
for k = 3:numel(fpth)
	addpath(genpath(fullfile(pth{1}, fpth{k})));
end
clear fpth k op


%% Set file names
loadfile = fullfile(pth{5},"all_demeaned_iteration2.mat");
fileName = "subjects_FBIRNmodes_demeaned";


%% Load and sort data

% Load formatted dFNC data
load(fullfile(pth{2}, "Data", "FBIRN_DFNC_table.mat"));
load(fullfile(pth{2}, "Data", "head_motion_meanFD.mat"));

% Confirm that IDs, data are properly indexed
assert(all(str2double(string(cell2mat(analysis_ID))) == str2double(analysis_data.Properties.RowNames)), "Data labels are not properly ordered!");
clear analysis_ID analysis_SCORE
assert(all(strcmpi(string(FILE_ID), string(analysis_data.Properties.VariableNames))), "Clinical variables are not properly ordered!");
clear FILE_ID

% add head motion to data array
analysis_data = [table(head_motion_meanFD, 'VariableNames',"Mean Head Motion"), analysis_data];
clear head_motion_meanFD


%% Convert table variables (should be placed in separate script)

% Replace numeric missing data code with NaN
analysis_data{:,:}(analysis_data{:,:} == -9999) = NaN;

% Set diagnosis, gender labels
labels.diagnosis = ["SZ"; "HC"];
labels.gender = ["M"; "F"];
labels.data = ["Diagnosis"; "Gender"];
labels.methods = ["Standard"; "Exact"];

% Identify table variables to change
i(1,:) = contains(analysis_data.Properties.VariableNames, "diagnosis");
i(2,:) = contains(analysis_data.Properties.VariableNames, "gender");

% generate string arrays
groups = labels.diagnosis(analysis_data{:,i(1,:)});
gender = labels.gender(analysis_data{:,i(2,:)});

% Convert variable type
analysis_data = convertvars(analysis_data, ["diagnosis(1:sz; 2:hc)","gender(1:male; 2:female)"], "string");

% Replace table numeric indices with strings
analysis_data{:,i(1,:)} = groups;
analysis_data{:,i(2,:)} = gender;
clear groups gender i

% Rename table variables
analysis_data = renamevars(analysis_data, ["age" "diagnosis(1:sz; 2:hc)" "gender(1:male; 2:female)"], ["Age" "Diagnosis" "Gender"]);


%% Set region labels & maps

% Load functional network labels
labels.FNC = readtable(fullfile(pth{4}, "NeuroMark_FNC_labels.xlsx")); % NeuroMark functional network labels & locations
labels.FNC = renamevars(labels.FNC, "SelectedComponentsAsRegionsOfInterest", "Functional Networks");

% Remove borders between functional domains
[r,~] = find(strcmpi(labels.FNC{:,:}, ""));   % find rows which separate functional domains
r = unique(r);
labels.ROI = array2table(cellfun(@str2num, labels.FNC{:,2:end}, 'UniformOutput', false), 'RowNames',labels.FNC{:,1}, 'VariableNames',labels.FNC.Properties.VariableNames(2:end));
labels.ROI(r,:) = [];

% Set functional domain labels
labels.FDs = labels.FNC(r,"Functional Networks");
labels.FDs = renamevars(labels.FDs, "Functional Networks", "Functional Domains");

% Set number of ROIs, FDs
N.ROI = size(labels.ROI,1);
N.FD = size(labels.FDs,1);

% Establish FN-level map of FDs
r = [r; size(labels.FNC,1)+1];
labels.FND = labels.FNC;
for i = 2:numel(r)
    labels.FND{r(i-1)+1:r(i)-1,"Functional Networks"} = repmat(labels.FDs{i-1,"Functional Domains"}, [r(i)-1-r(i-1) ,1]);
end
labels.FND(r(1:end-1),:) = [];

% Establish FN-level map of FDs
ind.FND = zeros(N.ROI, N.FD);
for d = 1:N.FD
    ind.FND(:,d) = strcmpi(labels.FND{:,"Functional Networks"}, labels.FDs{d,"Functional Domains"});
end
ind.FND = array2table(ind.FND, 'RowNames',labels.ROI.Properties.RowNames, 'VariableNames',labels.FDs{:,"Functional Domains"});
clear d i r


%% Concatenate and demean time series

% rename FNC data
FNC = cellfun(@transpose, DFNC_FBIRN, 'UniformOutput',false);

% demean FNC data
sFNC = cellfun(@mean, DFNC_FBIRN, 'UniformOutput',false);
sFNC = cellfun(@transpose, sFNC, 'UniformOutput',false);
FNC = cellfun(@minus, FNC, sFNC, 'UniformOutput',false);
sFNC = cell2mat(sFNC');


%% index time series

% Set counters
N.fig = 1;                                  % figures
N.TR = size(FNC{1},2);                      % time
N.conditions = numel(labels.diagnosis);     % conditions

% Index subjects
I = cell(1,N.conditions);
for j = 1:N.conditions
    I{j} = nnz(analysis_data{:,"Diagnosis"} == labels.diagnosis(j));
end
N.subjects = cell2table(I, 'VariableNames',labels.diagnosis);

% Locate & isolate site indices
ind.site = unique(analysis_data{:,'Site'});
clear I d j i

% convert variables to row form
I.subject = str2double(string(analysis_data.Properties.RowNames)');
I.diagnosis = analysis_data{:,"Diagnosis"}';
I.gender = analysis_data{:,"Gender"}';
I.site = analysis_data{:,"Site"}';
I.age = analysis_data{:,"Age"}';

% fill in data for each timepoint
f = fieldnames(I);
for j = 1:numel(f)
    I.(f{j}) = repmat(I.(f{j}), [N.TR 1]);
    I.(f{j}) = reshape(I.(f{j}), [sum(N.subjects{:,:})*N.TR 1]);    % reshape indices to column form
end
clear j f k

% Confirm that indices are in proper order
assert(all(unique(I.subject)==str2double(string(analysis_data.Properties.RowNames))), "Indices are out of order!");

% convert index to table
I = struct2table(I);


%% Set number of ICs

% % Find maximum number of components
% disp("Identifying number independent components from Marcenko-Pasteur distribution.");
% N.IC = NumberofIC(FNC.full);

% % Evaluate IC counts vs. captured variance
% d = factor(N.IC);
% N.IC = d(1)*d(2):d(1)*d(2):N.IC; clear d
% [ev, F(N.fig)] = evaluateICnumbers(N, ts);
% N.fig = N.fig + 1;

% Set number of modes to use in reconstruction
N.modes = 6;


%% Define filename based on parameters

% Get file list
fList = dir(fullfile(pth{5}, strcat(strjoin([fileName, "iteration"], '_'), '*.mat')));
fList = struct2table(fList);

% Set iteration number
a = false(size(fList,1),1);
for n = 1:size(fList,1)
    a(n) = matches("iteration", strsplit(string(fList{n,"name"}), '_'));
end
nIter = size(fList,1)-sum(a)+1;

% Set full filename
fileName = strjoin([fileName, strcat("iteration", num2str(nIter))], '_');
clear fList nIter a k n


%% Compute spectra

% Load dataset-level modes and time courses
load(loadfile, "Phi","mu","lambda","diagS","x0","dstnc");
Phi = Phi(:,:,strcmpi(labels.methods,"Exact"));

% Compute and plot mode spectra and mode power
[f, P, F(N.fig)] = DMD_spectrum(Phi, mu, 'plotit',1);   % power
F(N.fig).OuterPosition = [1 1 1055 1055]; hold on;      % increase figure size
title("DMD Power Spectrum");
xlim([min(f) max(f)]); ylim([0 max(P)]);

% sort power by frequency
[f_sort, i] = sort(f,1);
P_sort = P(i);

% Plot cumulative power of the modes
F(N.fig+1) = figure; F(N.fig+1).OuterPosition = [1 1 1055 1055];
plot(f_sort, cumsum(P_sort)./max(cumsum(P_sort),[],'all')); hold on;
plot([min(f_sort) max(f_sort)], [0.9 0.9], '-r');
xlabel("frequency (Hz)"); ylabel("% Cumulative Power");
xlim([min(f_sort) max(f_sort)]);
title("Cumulative DMD Power");
legend("Cumulative Power", "90% of Power", 'Location','southeast');

N.fig = N.fig+2;
clear i P_sort f_sort


%% Plot eigenvalues on unit circle

% Separate eigenvalues into real, imaginary parts
i = imag(lambda);
r = real(lambda);
c(:,1) = abs(lambda) > 1;
c(:,2) = abs(lambda) < 1;
c(:,3) = abs(lambda) == 1;

% Plot eigenvalues on unit circle
F(N.fig) = figure; N.fig = N.fig+1;
F(N.fig-1).OuterPosition = [1 1 1055 1055];
pbaspect([1 1 1]); hold on
s(1) = scatter(r(c(:,1)), i(c(:,1)), 'r',"filled");
s(2) = scatter(r(c(:,2)), i(c(:,2)), 'b',"filled");
s(3) = scatter(r(c(:,3)), i(c(:,3)), 'g',"filled");
% scatter(r, i, 'MarkerEdgeColor','k');

% Plot unit circle in real, imaginary space
theta = 0:0.1:2*pi+0.1;
x = cos(theta); y = sin(theta);
plot(x, y, '-k');
xlabel("Real"); ylabel("Imaginary");
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
legend(s, {'\lambda > 1', '\lambda < 1', '\lambda = 1'});
title("Dataset Eigenvalues");

clear i r c theta


%% Map group-level modes to subject time courses

% Preallocate arrays
tc = nan(N.modes, N.TR, sum(N.subjects{:,:}));
sm = nan(N.ROI*(N.ROI-1)/2, N.modes, sum(N.subjects{:,:}));

% Select modes to analyze
[~, ia, ~] = unique(P);
ia = flip(ia);  % list modes by power

% Regress subject-level modes and time courses from dataset-level modes
for s = 1:sum(N.subjects{:,:})
    for t = 1:N.TR
        tc(:,t,s) = regress(FNC{s}(:,t), Phi(:, ia(1:N.modes)));
    end
    for c = 1:N.ROI*(N.ROI-1)/2
        sm(c,:,s) = regress(FNC{s}(c,:)', tc(:,:,s)');
    end
end
clear m s t ans tc_i tc_r spatial_maps_real spatial_maps_imag

% % correct scaling of spatial maps and time courses
% spatial_maps = spatial_maps/10^2;
% tc = tc/10^2;


%% Visualize most powerful modes in random selection of subjects

% Select subjects to plot
sp = round(sum(N.subjects{:,:})*rand([2 1]));

% Plot real parts of most powerful modes for selected subjects
for c = 1:N.modes
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];

    subplot(2,3, 1);
    display_FNC(icatb_vec2mat(real(Phi(:,ia(c)))), [0.05 1.5]); hold on
    title(strjoin(["Mode ", num2str(ia(c)), ", Real Part"], ""));
    subplot(2,3, 4);
    display_FNC(icatb_vec2mat(imag(Phi(:,ia(c)))), [0.05 1.5]); hold on
    title(strjoin(["Mode ", num2str(ia(c)), ", Imaginary Part"], ""));

    for s = 1:numel(sp)
        subplot(2,3, s+1);
        display_FNC(icatb_vec2mat(real(sm(:,c,sp(s)))), [0.05 1.5]); hold on
        title(strjoin(["Subject ", num2str(sp(s)), ", Real Part"], ""));

        subplot(2,3, s+4);
        display_FNC(icatb_vec2mat(imag(sm(:,c,sp(s))')), [0.05 1.5]); hold on
        title(strjoin(["Subject ", num2str(sp(s)), ", Imaginary Part"], ""));

        sgtitle(strjoin(["Mode", num2str(ia(c))]));
    end
end
clear c s


%% Extract subject-level power spectra

P_sub = nan(N.modes, sum(N.subjects{:,:}));
for s = 1:sum(N.subjects{:,:})
    [~, P_sub(:,s)] = DMD_spectrum(sm(:,:,s), mu, 'plotit',0);   % power
end

% extract group means, medians, stds of power spectra
P_mean = nan(N.modes, N.conditions);
P_median = nan(N.modes, N.conditions);
P_std = nan(N.modes, N.conditions);
P_grp = nan(N.modes, max(N.subjects{:,:}), N.conditions);
for c = 1:N.conditions
    P_mean(:,c) = mean(P_sub(:, strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c))), 2);
    P_median(:,c) = median(P_sub(:, strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c))), 2);
    P_std(:,c) = std(P_sub(:, strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c))), 0, 2);
    P_grp(:, 1:N.subjects{:,labels.diagnosis(c)}, strcmpi(unique(analysis_data{:,"Diagnosis"}), ...
        labels.diagnosis(c))) = P_sub(:, strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c)));
end

% Select modes to analyze
[~, ia, ~] = unique(P);
ia = flip(ia);  % list modes by power

% plot group means and standard deviations of power spectra
F(N.fig) = figure; N.fig = N.fig+1;
F(N.fig-1).OuterPosition = [1 1 1920 1055];
subplot(3,1,1);
bar(f(ia(1:N.modes)), log(P_mean));
% errorbar(repmat(f(ia), [1 2]), P_mean(ia,:), P_std(ia,:), "_");
title("Means of Power Spectra"); legend(labels.diagnosis);
hold on;

% plot group medians and standard deviations of power spectra
subplot(3,1,2);
bar(f(ia(1:N.modes)), log(P_median));
% errorbar(repmat(f(ia), [1 2]), P_median(ia,:), P_std(ia,:), "_");
title("Medians of Power Spectra"); legend(labels.diagnosis);
hold on;

% plot group standard deviations of power spectra
subplot(3,1,3);
bar(f(ia(1:N.modes)), log(P_std));
title("Standard Deviations of Power Spectra"); legend(labels.diagnosis);
hold on;

% plot group power spectra distributions as violin plots
F(N.fig) = figure; N.fig = N.fig+1;
F(N.fig-1).OuterPosition = [1 1 1920 1055];
m = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
h = daviolinplot(P_sub', 'groups',analysis_data{:,"Diagnosis"}, 'xtlabels',f(ia(1:N.modes)), 'colors',m, 'boxcolors','same');   % each frequency is a different condition;
title("Power Spectra Distributions"); legend(labels.diagnosis);                                                                 % groups defined along subjects
hold on; clear m

% Run group comparisons between power spectra
p.ps = nan(N.modes,1);
for m = 1:N.modes
    [~, p.ps(m)] = kstest2(squeeze(P_grp(m,:,1)), squeeze(P_grp(m,:,2)));
end
h.ps = (p.ps < 0.05/N.modes);


%% Compare subject mode time series to frequency harmonic

% run periodogram analysis of each subject
pxx = nan(256, N.modes, sum(N.subjects{:,:}));
fc = nan(256, sum(N.subjects{:,:}));
for s = 1:sum(N.subjects{:,:})
    [pxx(:,:,s), fc(:,s)] = periodogram(squeeze(tc(:,:,s))');
end

% Plot real parts of most powerful modes for selected subjects
for c = 1:N.modes
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];
    for s = 1:numel(sp)
        subplot(numel(sp),1, s);
        scatter(1:N.TR, tc(c,:,sp(s))); hold on
        % plot(1:N.TR, max(tc(c,:,sp(s)),[],"all")*sin(real((1:N.TR)*2*f(ia(c)))));
        xlabel("Time Windows (2 sec)");
        title(strjoin(["Subject", num2str(sp(s))]));
        sgtitle([strjoin(["Mode", num2str(ia(c))]), strjoin([num2str(f(ia(c))), "Hz"])]);
    end
end
clear c s


%% view z-scored modes from random subjects

% organize spatial maps by group
sm_sep = cell(N.conditions,1);
i = nan(sum(N.subjects{:,:}), N.conditions);
sp = nan(2, N.conditions);
for c = 1:N.conditions
    sm_sep{c} = zscore(sm(:,:,strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c))));
    sm_sep{c} = permute(sm_sep{c}, [3 1 2]);

    i(:,c) = strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c));
    k = find(strcmpi(analysis_data{:,"Diagnosis"},labels.diagnosis(c)));
    sp(:,c) = k(round(rand(2,1)*N.subjects{:,labels.diagnosis(c)}));
end

for m = 1:N.modes
    for c = 1:N.conditions
        F(N.fig) = figure; N.fig = N.fig+1;
        F(N.fig-1).OuterPosition = [1 1 1920 1080];
        subplot(2,3, 1);
        display_FNC(icatb_vec2mat(real(squeeze(mean(sm_sep{c}(:,:,m),1)))'), [0.05 1.5]); hold on
        title(strjoin(["Mode", num2str(ia(m))]), "Real Part");
        subplot(2,3, 4);
        display_FNC(icatb_vec2mat(imag(squeeze(mean(sm_sep{c}(:,:,m),1)))'), [0.05 1.5]); hold on
        title(strjoin(["Mode", num2str(ia(m))]), "Imaginary Part");
    
        for s = 1:numel(sp(:,c))
            subplot(2,3, s+1);
            display_FNC(icatb_vec2mat(real(sm(:,m,sp(s,c)))'), [0.05 1.5]); hold on
            title(strjoin(["Subject ", num2str(sp(s,c)), ", Real Part"], ""));
    
            subplot(2,3, s+4);
            display_FNC(icatb_vec2mat(imag(sm(:,m,sp(s,c)))'), [0.05 1.5]); hold on
            title(strjoin(["Subject ", num2str(sp(s,c)), ", Imaginary Part"], ""));
        end
        sgtitle([strjoin(["z-Scored Mode", num2str(ia(m))]), labels.diagnosis(c)]);
    end
end

% display mean modes from each condition
for m = 1:N.modes
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1080 1080];
    for c = 1:N.conditions
        s = icatb_vec2mat(squeeze(mean(sm_sep{c},1))');
        subplot(2,2,c);
        display_FNC(real(squeeze(s(m,:,:))), [0.05 1.5]); hold on
        title(labels.diagnosis(c), "Real Part");
        subplot(2,2,c+2);
        display_FNC(imag(squeeze(s(m,:,:))), [0.05 1.5]); hold on
        title(labels.diagnosis(c), "Imaginary Part")
    end
    sgtitle(["Mean z-Scored FNC", strjoin(["Mode", num2str(ia(m))])]);
end
clear c s m k


%% Search for group-level changes in spatial maps

% two-sample t-test
[~,p.t2] = ttest2(sm_sep{1}, sm_sep{2});
p.t2 = squeeze(p.t2);

% Kolmogorov-Smirnov test
for m = 1:N.modes
    for n = 1:(N.ROI*(N.ROI-1)/2)
        [~,p.ks(n,m,1)] = kstest2(real(squeeze(sm_sep{1}(:,n,m))), real(squeeze(sm_sep{2}(:,n,m))));
        [~,p.ks(n,m,2)] = kstest2(imag(squeeze(sm_sep{1}(:,n,m))), imag(squeeze(sm_sep{2}(:,n,m))));
    end
end

% Bonferroni correction
h.BF = (p.ks < 0.05/(N.ROI*(N.ROI-1)/2));
h.BF = h.BF(:,:,1) | h.BF(:,:,2);

% Benjamini-Hochberg FDR
h.BH = zeros(N.ROI*(N.ROI-1)/2, N.modes, 2);
for c = 1:2
    hd = zeros(N.ROI*(N.ROI-1)/2, N.modes);
    pval = reshape(p.ks(:,:,c), [numel(p.ks(:,:,c)) 1]);
    rejectedH0s = FDR_benjHoch(pval, 0.05, 'positive', true, false, false);
    hd(rejectedH0s) = 1;
    h.BH(:,:,c) = hd;
end
h.BH = h.BH(:,:,1) | h.BH(:,:,2);

% MA FDR
h.MA = zeros(N.ROI*(N.ROI-1)/2, N.modes, 2);
for c = 1:2
    hd = zeros(N.ROI*(N.ROI-1)/2, N.modes);
    pval = reshape(p.ks(:,:,c), [numel(p.ks(:,:,c)) 1]);
    fdr = mafdr(pval);
    hd(fdr<0.05) = 1;
    h.MA(:,:,c) = hd;
end
h.MA = h.MA(:,:,1) | h.MA(:,:,2);
clear c m n hd rejectedH0s pval


%% Visualize significant connections

% Bonferroni-significant connections
if nnz(h.BF) > 0
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];
    for m = 1:N.modes
        if nnz(h.BF(:,m)) > 0
            subplot(2,3,m);
            f = icatb_vec2mat(h.BF(:,m));
            [r,c] = find(triu(f));
            sm_mask = tril(icatb_vec2mat(zscore(Phi(:,ia(m)))));
            display_FNC(real(sm_mask), [0.05 1.5]); hold on
            scatter(c, r, 30, 'r', "square", "filled"); hold on
            title(strjoin(["Significant Connections of Mode", num2str(ia(m))]));
        end
    end
    sgtitle("Bonferroni Correction");
end
clear f sm_mask s m c r

% Benjamini-Hochberg significant connections
if nnz(h.BH) > 0
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];
    for m = 1:N.modes
        if nnz(h.BH(:,m)) > 0
            subplot(2,3,m);
            subplot(2,3,m);
            f = icatb_vec2mat(h.BH(:,m));
            [r,c] = find(triu(f));
            sm_mask = tril(icatb_vec2mat(zscore(Phi(:,ia(m)))));
            display_FNC(real(sm_mask), [0.05 1.5]); hold on
            scatter(c, r, 30, 'r', "square", "filled"); hold on
            title(strjoin(["Significant Connections of Mode", num2str(ia(m))]));
        end
    end
    sgtitle("Benjamini-Hochberg Correction");
end
clear f sm_mask s m c r

% MA-FDR significant connections
if nnz(h.MA) > 0
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];
    for m = 1:N.modes
        if nnz(h.MA(:,m)) > 0
            subplot(2,3,m);
            subplot(2,3,m);
            f = icatb_vec2mat(h.MA(:,m));
            [r,c] = find(triu(f));
            sm_mask = tril(icatb_vec2mat(zscore(Phi(:,ia(m)))));
            display_FNC(real(sm_mask), [0.05 1.5]); hold on
            scatter(c, r, 30, 'r', "square", "filled"); hold on
            title(strjoin(["Significant Connections of Mode", num2str(ia(m))]));
        end
    end
    sgtitle("False Discovery Rate Correction");
end
clear f sm_mask s m c r


%% Network-Based Statistic

% Calculate the NBS
nbs = cell(N.modes, 1);
STATS = cell(N.modes, 1);
GLM = cell(N.modes, 1);
storarray = cell(N.modes, 1);
% i = horzcat(ones(sum(N.subjects{:,:}), 1), i);
i = array2table(i, "VariableNames", labels.diagnosis);
contrast = [1 -1; -1 1];
tstat = [3 3.5];
for m = 1:N.modes
    mlabel(m) = strjoin(["Mode", num2str(ia(m))]);
    EC = squeeze(sm(:,m,:));
    EC = permute(icatb_vec2mat(EC'), [2 3 1]);
    [nbs{m}, STATS{m}, GLM{m}, storarray{m}] = runNBS(EC, contrast, i, N, tstat, labels);
end
nbs = cell2table(nbs', "VariableNames",mlabel);
clear c m s mlabel

% Display the NBS
col = ["r" "b"];
l = find(~cellfun(@isempty, storarray));
for t = 1:numel(tstat)
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1920 1080];
    for m = 1:length(l)
        subplot(2,3,l(m)); pbaspect([1 1 1]);
        display_FNC(zeros(N.ROI), [0.05 1.5], [], false); hold on
        for c = 1:2*nchoosek(2,2)
            [r, cl] = find(full(nbs.(l(m)){1,1}.(c){t,1}));
            k(c) = scatter(cl, r, 30, col(c), "square", "filled"); hold on
        end
        legend(k, nbs.(l(m)){1,1}.Properties.VariableNames);
        title(strjoin(["Significant Connections of Mode", num2str(ia(l(m)))]));
    end
    sgtitle(["Network-Based Statistic", strjoin(["t-statistic:", num2str(tstat(t))])]);
end
clear m c col r cl k l t


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(pth{5}, fileName), 'compact');
for c = 1:N.fig-1
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'svg');
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'jpeg');
end
clear c F a ax axes ts i ans

% Save files
N.fig = N.fig - 1;
save(fullfile(pth{5}, fileName));
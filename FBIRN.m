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

% Add relevant paths
addpath(genpath(fullfile(pth{2}, fpth{1})));
addpath(fullfile(pth{1}, fpth{2}));
for k = 3:numel(fpth)
	addpath(genpath(fullfile(pth{1}, fpth{k})));
end
clear fpth k op


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

% % Sort data by site
% [~, I(:,1)] = sort(analysis_data{:,"Site"});
% analysis_data = analysis_data(I,:);
% DFNC_FBIRN = DFNC_FBIRN(I,:);
% clear I

% Set diagnosis, gender labels
labels.diagnosis = ["SZ";"HC"];
labels.gender = ["M";"F"];
labels.data = ["Diagnosis", "Gender"];


%% Index counters

% Set figure counter
N.fig = 1;

% Index time
N.TR = size(DFNC_FBIRN{1},1);

% identify diagnosis column
i = contains(analysis_data.Properties.VariableNames, 'diagnosis');

% Index conditions
d = unique(analysis_data{:,i});
N.conditions = numel(d);

% Index subjects
I = cell(1,N.conditions);
for j = 1:N.conditions
    I{j} = nnz(analysis_data{:,i} == d(j));
end
N.subjects = cell2table(I, 'VariableNames',labels.diagnosis);

% Locate & isolate site indices
ind.site = unique(analysis_data{:,'Site'});
clear I d j i


%% Convert table variables (should be placed in separate script)

% Replace numeric missing data code with NaN
analysis_data{:,:}(analysis_data{:,:} == -9999) = NaN;

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


%% Concatenate and index time series

% rename FNC data
FNC.subj = cellfun(@transpose, DFNC_FBIRN, 'UniformOutput',false);
FNC.full = cell2mat(DFNC_FBIRN)';
clear DFNC_FBIRN

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

% Set number of ICs
% N.IC = 8;


%% Define filename based on parameters

% Get number of ICs
fileName = "test"; % strjoin(["DMD", strcat(num2str(N.IC), "ICs")], "_");

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


%% Isolate subject-level components & activity from dFNC

% Preallocate arrays
Phi.subj = nan(N.ROI*(N.ROI-1)/2, N.TR-1, sum(N.subjects{:,:}), 2);
mu.subj = nan(N.TR-1, sum(N.subjects{:,:}));
lambda.subj = nan(N.TR-1, sum(N.subjects{:,:}));
diagS.subj = nan(N.TR-1, sum(N.subjects{:,:}));
x0.subj = nan(N.ROI*(N.ROI-1)/2, sum(N.subjects{:,:}));
dstnc.subj = cell(sum(N.subjects{:,:}), 1);

% Run subject-level DMD
for s = 1:sum(N.subjects{:,:})
    % Generate subject-level X, Y matrices
    X = FNC.subj{s}(:, 1:N.TR-1);
    Y = FNC.subj{s}(:, 2:N.TR);

    % Run DMD
    [Phi.subj(:,:,s,1), mu.subj(:,s), lambda.subj(:,s), diagS.subj(:,s), x0.subj(:,s)] = DMD(X, Y, 'dt',2);  % standard
    [Phi.subj(:,:,s,2), ~, ~, ~, ~] = DMD(X, Y, 'dt',2, 'exact',true);                   % exact

    % Compare exact vs. standard DMD
    D = Phi.subj(:,:,s,1) - Phi.subj(:,:,s,2);
    d = nnz(abs(D) <= eps);
    if ~d
        dstnc.subj{s} = [];
    else
        dstnc.subj{2} = abs(Phi.subj(:,:,2,1) - Phi.subj(:,:,2,2));
    end
end

% Check if exact and standard SVD produce same outputs
i = cellfun(@isempty, dstnc.subj);
if nnz(i) == length(dstnc.subj)
    Phi.subj = squeeze(Phi.subj(:,:,:,2));    % keep exact DMD
    dstnc = rmfield(dstnc, "subj");
end
clear s X Y d D i


%% Isolate group-level components & activity from dFNC

% Preallocate arrays
Phi.grp = nan(N.ROI*(N.ROI-1)/2, N.ROI*(N.ROI-1)/2, N.conditions, 2);
mu.grp = nan(N.ROI*(N.ROI-1)/2, N.conditions);
lambda.grp = nan(N.ROI*(N.ROI-1)/2, N.conditions);
diagS.grp = nan(N.ROI*(N.ROI-1)/2, N.conditions);
x0.grp = nan(N.ROI*(N.ROI-1)/2, N.conditions);
dstnc.grp = cell(N.conditions, 1);

% Run group-level DMD
for g = 1:N.conditions
    % Generate group-level X, Y matrices
    X = FNC.subj(analysis_data{:,'Diagnosis'} == labels.diagnosis(g))';
    Y = FNC.subj(analysis_data{:,'Diagnosis'} == labels.diagnosis(g))';
    for s = 1:numel(X)
        X{s} = X{s}(:, 1:N.TR-1);
        Y{s} = Y{s}(:, 2:N.TR);
    end
    X = cell2mat(X);
    Y = cell2mat(Y);

    % Run DMD
    [Phi.grp(:,:,g,1), mu.grp(:,g), lambda.grp(:,g), diagS.grp(:,g), x0.grp(:,g)] = DMD(X, Y, 'dt',2);  % standard
    [Phi.grp(:,:,g,2), ~, ~, ~, ~] = DMD(X, Y, 'dt',2, 'exact',true);                   % exact

    % Compare exact vs. standard DMD
    D = Phi.grp(:,:,g,1) - Phi.grp(:,:,g,2);
    d = nnz(abs(D) <= eps);
    if ~d
        dstnc.grp{g} = [];
    else
        dstnc.grp{g} = abs(Phi.grp(:,:,g,1) - Phi.grp(:,:,g,2));
    end
end

% Check if exact and standard SVD produce same outputs
i = cellfun(@isempty, dstnc.grp);
if nnz(i) == length(dstnc.grp)
    Phi.grp = squeeze(Phi.grp(:,:,:,2));    % keep exact DMD
    dstnc = rmfield(dstnc, "grp");
end
clear g X Y d D i


%% Compute spectra

% compute DMD spectrum for each subject
[f, P, F] = DMD_spectrum(Phi(:,:,1), mu(:,1), 'plotit',1);  % power
F(numel(F)).OuterPosition = [1 1 1440 1055]; hold on;   % increase figure size
plot(sort(f), 1./(sort(f)), '--g');
title("Mode Power Spectrum"); ylim([0 max(P)]);
legend({'Frequency Power', '$\frac{1}{f}$'}, 'Interpreter','latex');

% Compute phases
phi = atan(imag(Phi(:,:,1))./real(Phi(:,:,1)));

% Compute cumulative power of the modes
[f_sort, i] = sort(f,1);
P_sort = P(i);
lambda_sort = lambda(i,:);
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1055 1055];
plot(f_sort, cumsum(P_sort)./max(cumsum(P_sort),[],'all')); hold on;
plot([min(f_sort) max(f_sort)], [0.9 0.9], '-r');
xlabel("frequency (Hz)"); ylabel("% Cumulative Power");
title("Cumulative Power of the Frequencies");
legend("Cumulative Sum", "90% of Power");
clear P_sort f_sort


%% 

% Compare subject modes to group modes (confirm in same order)
dstnc.svg = cell(sum(N.subjects{:,:}), 1);
for s = 1:sum(N.subjects{:,:})
    i = find(strcmpi(analysis_data{7,'Diagnosis'}, labels.diagnosis));
    dstnc.svg{s} = abs(Phi.subj(:,:,s) - Phi.grp(:,1:N.TR-1,i));
end
clear i


%% Check reconstruction error

% Set number of modes to use in reconstruction
N.modes = 3;

% Compare method reconstruction
Phi_sort = Phi(:,i,1);

% Open figure (large)
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];

% Compute MSE per sample (partial)
[Xhat.partial, ~] = DMD_recon(Phi_sort(:,1:N.modes,1), lambda_sort(1:N.modes,1), x0(:,1), N.TR);    % first five modes
e = FNC.subj{1}-Xhat.partial;
msqe.partial = abs(sum(e.^2))/(N.ROI*(N.ROI-1)/2);

% Visualize estimated sFNC (partial)
subplot(2,3,1);
display_FNC(real(icatb_vec2mat(mean(Xhat.partial,2))), [0.25 1.5]);
title(strjoin(["Reconstructed sFNC (first", num2str(N.modes), "modes)"])); hold on;

% Visualize MSE per sample (partial reconstruction)
subplot(2,3,4);
stem(msqe.partial); axis tight; hold on
xlabel('samples'); ylabel('MSE');
title(strjoin(["MSE per sample (first", num2str(N.modes), "modes)"]));

% Compute MSE per sample (full)
[Xhat.full, ~] = DMD_recon(Phi_sort(:,:,1), lambda_sort(:,1), x0(:,1), N.TR);
e = FNC.subj{1}-Xhat.full;
msqe.full = abs(sum(e.^2))/(N.ROI*(N.ROI-1)/2);

% Visualize estimated sFNC (full)
subplot(2,3,2);
display_FNC(real(icatb_vec2mat(mean(Xhat.full,2))), [0.25 1.5]);
title("Reconstructed sFNC (full)"); hold on;

% Visualize MSE per sample (full reconstruction)
subplot(2,3,5);
stem(msqe.full); axis tight; hold on
xlabel('samples'); ylabel('MSE');
title("MSE per sample (full reconstruction)");

% Visualize actual sFNC
subplot(2,3,6);
display_FNC(icatb_vec2mat(mean(FNC.subj{1},2)), [0.25 1.5]); hold on;
title("Original sFNC"); hold on;

% Display difference between partial and full reconstructions
subplot(2,3,3);
display_FNC(real(icatb_vec2mat(mean(Xhat.full,2) - mean(Xhat.partial,2))), [0.25 1.5]); hold on;
title('Full - Partial Reconstruction');


%% Plot eigenvalues on unit circle

% Separate eigenvalues into real, imaginary parts
i = imag(lambda(:,1));
r = real(lambda(:,1));
c(:,1) = abs(lambda(:,1)) > 1;
c(:,2) = abs(lambda(:,1)) < 1;
c(:,3) = abs(lambda(:,1)) == 1;

% Plot eigenvalues on unit circle
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1055 1055]; hold on
s(1) = scatter(r(c(:,1)), i(c(:,1)), 'MarkerFaceColor','r');
s(2) = scatter(r(c(:,2)), i(c(:,2)), 'MarkerFaceColor','b');
s(3) = scatter(r(c(:,3)), i(c(:,3)), 'MarkerFaceColor','g');
% scatter(r, i, 'MarkerEdgeColor','k');

% Plot unit circle in real, imaginary space
theta = 0:0.1:2*pi+0.1;
x = cos(theta); y = sin(theta);
plot(x, y, '-k');
xlabel("Real"); ylabel("Imaginary");
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
legend(s, {'\lambda > 1', '\lambda < 1', '\lambda = 1'});


%% Plot FN time courses from single module (per mode)

% sort frequencies
[f_sort, i] = unique(f);
Phi_sort = Phi(:,i,1);
lambda_sort = lambda(i,1);

% set mask
r = [6 6 6 7 7 7];
c = [2 3 4 2 3 4];
ind.rc = horzcat(r', c');

% convert masks to linear indices
m = zeros(N.ROI, N.ROI);
m(r,c) = 1;
ind.lin = find(m);

% Plot mask
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
subplot(2,4,1);
imagesc(m); colormap bone; colorbar; pbaspect([1 1 1]);
title("Timecourse Mask");
xlabel("Neuromark Functional Networks");
ylabel("Neuromark Functional Networks");

% Plot original FN courses over time
subplot(2,4,5);
plot(1:N.TR, FNC.subj{1}(ind.lin(:,1),:)); hold on;
title(strjoin("Original FNC Values"));
xlabel("Time Points"); ylabel("Real Amplitude");
legend(num2str(ind.rc));

% Plot FN courses over time as a function of number of modes
ii = [2 3 4 6 7 8];
for k = 1:length(c)
    [Xhat, ~] = DMD_recon(Phi_sort(:,i(k+1),1), lambda_sort(i(k+1),1), x0(:,1), N.TR);    % first five modes
    subplot(2,4,ii(k));
    plot(1:N.TR, Xhat(ind.lin,:)); hold on;
    title("Reconstructed FNC Values", strjoin(["f =", num2str(f_sort(k+1)), "Hz"]));
    xlabel("Time Points"); ylabel("Real Amplitude");
    legend(num2str(ind.rc));
end
clear i ii ind r c k m n e Xhat lambda_sort


%% Visualize modes

% remove duplicate (negative) modes
[f_sort, i, irev] = unique(f);
P_sort = P(i);
Phi_sort = Phi(:,i,1);
phi_sort = phi(:,i,1);

% get amplitudes as function of frequency
l.r = max(abs(real(Phi_sort(:,:,1))));
l.i = max(abs(imag(Phi_sort(:,:,1))));

% Visualize static mode
Phi_mat = icatb_vec2mat(squeeze(Phi_sort(:,1,1)));
phase_mat = icatb_vec2mat(squeeze(phi_sort(:,1,1)));
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
subplot(1,2,1); display_FNC(real(Phi_mat), [0.25 1.5]); title("Mode (Real Part)"); hold on;
subplot(1,2,2); display_FNC(imag(Phi_mat), [0.25 1.5], [-max(abs(l.i)) max(abs(l.i))]); title("Mode (Imaginary Part)"); hold on;
sgtitle(strjoin(["f =" , num2str(f_sort(1))]));
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
display_FNC(real(phase_mat), [0.25 1.5], [-pi/2 pi/2]); title('Phases at \omega = 0'); hold on;

% visualize dominant harmonic modes
for j = 2:nnz(cumsum(P_sort)./max(cumsum(P_sort)) < 0.9)
    Phi_mat = icatb_vec2mat(squeeze(Phi_sort(:,j,1)));
    phase_mat = icatb_vec2mat(squeeze(phi_sort(:,j,1)));
    F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
    subplot(1,2,1); display_FNC(real(Phi_mat), [0.25 1.5]); title("Mode (Real Part)"); hold on;
    subplot(1,2,2); display_FNC(imag(Phi_mat), [0.25 1.5]); title("Mode (Imaginary Part)"); hold on;
    sgtitle(strjoin(["f =" , num2str(f_sort(j))]));

    F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
    display_FNC(real(phase_mat), [0.25 1.5], [-pi/2 pi/2]); title(strjoin({'Phases at \omega = ', num2str(2*f_sort(j)), '\pi'}, '')); hold on;
end

% visualize amplitudes as function of frequency
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
plot(f_sort, l.r, 'r'); hold on
plot(f_sort, l.i, 'b');
title('Absolute Amplitudes by Frequency');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
legend('Real','Imaginary');

clear i j Phi_mat phase_mat Phi_sort phi_sort f_sort P_sort


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(pth{5}, fileName), 'compact');
for c = 1:numel(F)
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'svg');
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'jpeg');
end
clear c F a ax axes ts

% Save files
N.fig = N.fig - 1;
save(fullfile(pth{5}, fileName));
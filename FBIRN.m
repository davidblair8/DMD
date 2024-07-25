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
labels.FNC = readtable(fullfile(pth{2}, "Data", "NeuroMark_FNC_labels.xlsx")); % NeuroMark functional network labels & locations
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
fList = dir(fullfile(pth{4}, strcat(strjoin([fileName, "iteration"], '_'), '*.mat')));
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


%% Isolate components & activity from dFNC

% extract components & activities from dFNC
SVD = nan(N.ROI*(N.ROI-1)/2, N.TR-1, sum(N.subjects{:,:}));
Hu = nan(N.ROI*(N.ROI-1)/2, N.TR-1, sum(N.subjects{:,:}));
corrphi = nan(2, 2, sum(N.subjects{:,:}));
mu = nan(N.TR-1, sum(N.subjects{:,:}));
lambda = nan(N.TR-1, sum(N.subjects{:,:}));
diagS = nan(N.TR-1, sum(N.subjects{:,:}));
x0 = nan(N.ROI*(N.ROI-1)/2, sum(N.subjects{:,:}));
for s = 1:sum(N.subjects{:,:})
    [SVD(:,:,s), mu(:,s), lambda(:,s), diagS(:,s), x0(:,s)] = DMD(FNC.subj{s}, 'dt',1);
    [Hu(:,:,s), ~, ~, ~, ~] = DMD_Hu(FNC.subj{s}, 'dt',1);
    corrphi(:,:,s) = corrcoef(SVD(:,:,s), Hu(:,:,s));
end
clear s

% locate modal non-matches
f = nnz((abs(corrphi-1) > eps));
if f == 0
    Phi = Hu;
    clear SVD Hu
else
    error("DMD methods in conflict!");
end


%% Compute spectra

% compute eigenvalue spectra
spect.mu = mu.*conj(mu);                % the fourier spectrum of modes (mu = log(lambda)/dt)
spect.lambda = lambda.*conj(lambda);    % DMD spectrum of modes 
spect.Phi = Phi.*conj(Phi);             % the modes

% compute DMD spectrum
[f, P, F] = DMD_spectrum(Phi(:,:,1), mu(:,1), 'plotit',1);  % power
F(numel(F)).OuterPosition = [1 1 1440 1055]; hold on;   % increase figure size
title("Mode Power Spectrum");

% Compute phases
phi = atan(imag(Phi(:,:,1))./real(Phi(:,:,1)));

% Compute cumulative power of the modes
[fc, i] = sort(f,1);
Pc = P(i);
Phi_sort = Phi(:,i,1);
lambda_sort = lambda(i,:);
phi_sort = phi(:,i,1);
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1055 1055];
plot(fc, cumsum(Pc)./max(cumsum(Pc),[],'all')); hold on;
plot([min(fc) max(fc)], [0.9 0.9], '-r');
xlabel("frequency (Hz)"); ylabel("% Cumulative Power");
title("Cumulative Power of the Frequencie");
legend("Cumulative Sum", "90% of Power");
clear Pc fc


%% Check reconstruction error

% Set number of modes to use in reconstruction
N.modes = 3;

% Visualize actual sFNC
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1055 1055];
display_FNC(icatb_vec2mat(mean(FNC.subj{1},2)), [0.25 1.5]);
title("Original sFNC"); hold on;

% Compare method reconstruction
Phi_sort = Phi(:,i,1);

% Compute MSE per sample (full)
[Xhat, ~] = DMD_recon(Phi_sort(:,:,1), lambda_sort(:,1), x0(:,1), N.TR);
e = FNC.subj{1}-real(Xhat);
msqe.full = abs(sum(e.^2))/(N.ROI*(N.ROI-1)/2);

% Visualize estimated sFNC (full)
F(numel(F)+ 1) = figure;
F(numel(F)).OuterPosition = [1 1 1440 1055]; hold on;   % increase figure size
subplot(2,2,1);
display_FNC(icatb_vec2mat(real(mean(Xhat,2))), [0.25 1.5]);
title("Reconstructed sFNC (full)"); hold on;

% Visualize MSE per sample (full reconstruction)
subplot(2,2,3);
stem(msqe.full); axis tight; hold on
xlabel('samples'); ylabel('MSE');
title("Reconstructed sFNC (full)"); hold on;

% Compute MSE per sample (partial)
[Xhat, ~] = DMD_recon(Phi_sort(:,1:N.modes,1), lambda_sort(1:N.modes,1), x0(:,1), N.TR);    % first five modes
e = FNC.subj{1}-real(Xhat);
msqe.partial = abs(sum(e.^2))/(N.ROI*(N.ROI-1)/2);

% Visualize estimated sFNC (partial)
subplot(2,2,2);
display_FNC(icatb_vec2mat(real(mean(Xhat,2))), [0.25 1.5]);
title(strjoin(["Reconstructed sFNC (first", num2str(N.modes), "modes)"])); hold on;

% Visualize MSE per sample (partial reconstruction)
subplot(2,2,4);
stem(msqe.partial); axis tight; hold on
xlabel('samples'); ylabel('MSE');
title(strjoin(["MSE per sample (first", num2str(N.modes), "modes)"]));


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
xlim([-1 1]); ylim([-1 1]);
legend(s, {'\lambda > 1', '\lambda < 1', '\lambda = 1'});


%% Plot FN time courses from single module (per mode)

% set mask for mode 1 (omega = 0)
clear r c
r(1,:) = [4, 5, 6, 7];
c(1,:) = [3, 3, 3, 3];

% set mask for mode 2 (omega = 0.012811)
r(2,:) = [25, 25, 25, 25];
c(2,:) = [20, 21, 22, 23];

% set mask for mode 3 (omega = 0.024677)
r(3,:) = [8, 8, 8, 8];
c(3,:) = [3, 4, 5, 6];

% set mask for mode 4 (omega = 0.03742)
r(4,:) = [9, 9, 9, 9];
c(4,:) = [3, 4, 5, 6];

% set mask for mode 5 (omega = 0.052477)
r(5,:) = [11, 11, 11, 11];
c(5,:) = [3, 4, 5, 6];

% set mask for mode 6 (omega = 0.065494)
r(6,:) = [36, 36, 36, 36];
c(6,:) = [27, 28, 29, 30];

% convert masks to linear indices
m = zeros(size(c,1), N.ROI, N.ROI);
i = nan(size(c))';
for k = 1:size(c,1)
    m(k, r(k,:), c(k,:)) = 1;
    n = squeeze(m(k,:,:));
    e = tril(ones(N.ROI), -1);
    n = n(e==1);
    i(:,k) = find(n);
end

% Plot original FN courses over time
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1150 1055];
plot(1:N.TR, FNC.subj{1}(i(:,1),:)); hold on;
title(strjoin("Original FNC Values"));
xlabel("Time Points"); ylabel("Real Amplitude");
legend(num2str(i));

% Plot FN courses over time as a function of number of modes
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
for k = 1:size(c,1)
    [Xhat, ~] = DMD_recon(Phi_sort(:,k,1), lambda_sort(k,1), x0(:,1), N.TR);    % first five modes
    subplot(2,3,k);
    plot(1:N.TR, Xhat(i(:,k),:)); hold on;
    title(strjoin(["Reconstructed FNC Values for FNs at Mode", num2str(k)]));
    xlabel("Time Points"); ylabel("Real Amplitude");
    legend(num2str(i(:,k)));
end
clear i r c k m n e


%% Visualize modes

% remove duplicate (negative) modes
[f, i, irev] = unique(f);
P = P(i);
Phi_sort = Phi_sort(:,i,1);
phi_sort = phi_sort(:,i,1);

% get amplitudes as function of frequency
l.r = max(abs(real(Phi_sort(:,:,1))));
l.i = max(abs(imag(Phi_sort(:,:,1))));

% Visualize static mode
Phi_mat = icatb_vec2mat(squeeze(Phi_sort(:,1,1)));
phase_mat = icatb_vec2mat(squeeze(phi_sort(:,1,1)));
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
subplot(1,2,1); display_FNC(real(Phi_mat), [0.25 1.5]); title('Mode (Real Part) at \omega = 0'); hold on;
subplot(1,2,2); display_FNC(imag(Phi_mat), [0.25 1.5], [-max(abs(l.i)) max(abs(l.i))]); title('Mode (Imaginary Part) at \omega = 0'); hold on;
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
display_FNC(real(phase_mat), [0.25 1.5], [-pi/2 pi/2]); title('Phases at \omega = 0'); hold on;

% visualize dominant harmonic modes
for j = 2:nnz(cumsum(P)./max(cumsum(P)) < 0.9)
    Phi_mat = icatb_vec2mat(squeeze(Phi_sort(:,j,1)));
    phase_mat = icatb_vec2mat(squeeze(phi_sort(:,j,1)));
    F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1920 1055];
    subplot(1,2,1); display_FNC(real(Phi_mat), [0.25 1.5]); title(strjoin({'Mode (Real Part) at \omega = ', num2str(2*f(j)), '\pi'}, '')); hold on;
    subplot(1,2,2); display_FNC(imag(Phi_mat), [0.25 1.5]); title(strjoin({'Mode (Imaginary Part) at \omega = ', num2str(2*f(j)), '\pi'}, '')); hold on;

    F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
    display_FNC(real(phase_mat), [0.25 1.5], [-pi/2 pi/2]); title(strjoin({'Phases at \omega = ', num2str(2*f(j)), '\pi'}, '')); hold on;
end
clear i j Phi_mat phase_mat Phi_sort phi_sort

% visualize amplitudes as function of frequency
F(numel(F)+1) = figure; F(numel(F)).OuterPosition = [1 1 1100 1055];
plot(f, l.r, 'r'); hold on
plot(f, l.i, 'b');
title('Absolute Amplitudes by Frequency');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
legend('Real','Imaginary');

%% Visualize reconstruction using individual modes


% save results
savefig(F, fileName, 'compact');
clear F; % close all
save(fileName);


%% Test for group-level changes in power spectra


%% Regress power spectra against clinical variables


%% Visualise static FNCs and spectral power arrays

% import NeuroMark templates from GIFT toolbox
fname = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'Neuromark_fMRI_1.0.nii');

% convert network domain label format for connectogram
network_names = cell(numel(labels.FDs), 2);
for n = 1:numel(labels.FDs)
    network_names{n,1} = string(labels.FDs{n,:});
    network_names{n,2} = find(strcmpi(string(labels.FND{:,"Functional Networks"}), string(labels.FDs{n,:})))';
end

% Compile spectral power arrays



%% Visualize modes with significant group-level power alterations

% import NeuroMark templates from GIFT toolbox
fname = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'Neuromark_fMRI_1.0.nii');

% convert network domain label format for connectogram
network_names = cell(numel(labels.FDs), 2);
for n = 1:numel(labels.FDs)
    network_names{n,1} = string(labels.FDs{n,:});
    network_names{n,2} = find(strcmpi(string(labels.FND{:,"Functional Networks"}), string(labels.FDs{n,:})))';
end
clear n fname


%% Plot modes as combined plot

% get screen dimensions
dim = get(0,'screensize');

% Horizontal
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).OuterPosition = [floor(dim(3)/2) 0 floor(dim(3)/2) dim(4)];
F(N.fig-1).Color = 'black';
for j = 1:N.IC
    % Plot mode correlation matrices

    % Plot spectral power as violin plots

    % Plot mode connectograms
    
end

% Save figure


% Vertical
F(N.fig) = figure; N.fig = N.fig + 1;
F(N.fig-1).OuterPosition = get(0,'screensize'); F(N.fig-1).Color = 'black';
for j = 1:N.IC
    % Plot mode correlation matrices

    % Plot spectral power as violin plots

    % Plot mode connectograms
    
end

% Save figure

clear dim j

% Save figure
saveas(F(N.fig-1), fullfile("/Users/David/Library/CloudStorage/GoogleDrive-dblair@gsu.edu/My Drive/Calhoun/Results/Regression/Figures", strjoin([fileName, "vertical"], '-')), 'svg');
% saveas(F(N.fig-1), fullfile(pth{4}, "Figures", strjoin([fileName, "vertical"], '-')), 'svg');


%% Plot large mode matrices, connectograms, and power spectra

% get screen dimensions
dim = get(0,'screensize');

for j = 1:N.IC
    % source matrices
    

    % Plot entropies
    

    % Connectogram
    
end
clear k ts fname network_names n j c y dim


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(pth{4}, fileName), 'compact');
for c = 1:numel(F)
    saveas(F(c), fullfile(pth{3}, "Images", strjoin([fileName, num2str(c)], '-')), 'svg');
end
clear c F a ax axes

% Save files
N.fig = N.fig - 1;
clear ts;
save(fullfile(pth{4}, fileName));
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


%% Concatenate and index time series

% rename FNC data
FNC = cellfun(@transpose, DFNC_FBIRN, 'UniformOutput',false);

% demean FNC data
mFNC = cellfun(@mean, cellfun(@mean, DFNC_FBIRN, 'UniformOutput',false), 'UniformOutput',false);
mFNC = cellfun(@transpose, mFNC, 'UniformOutput',false);
FNC = cellfun(@minus, FNC, mFNC, 'UniformOutput',false);
clear DFNC_FBIRN mFNC

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
N.modes = 3;


%% Define filename based on parameters

% define core file name
fileName = "subjects_demeaned";

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


%% Isolate group-level components & activity from dFNC

% Preallocate arrays
Phi = nan(N.ROI*(N.ROI-1)/2, N.TR-1, sum(N.subjects{:,:}), numel(labels.methods));
mu = nan(N.TR-1, sum(N.subjects{:,:}));
lambda = nan(N.TR-1, sum(N.subjects{:,:}));
diagS = nan(N.TR-1, sum(N.subjects{:,:}));
x0 = nan(N.ROI*(N.ROI-1)/2, sum(N.subjects{:,:}));
dstnc = zeros(N.ROI*(N.ROI-1)/2, N.TR-1, sum(N.subjects{:,:}));
sFNC = nan(N.ROI*(N.ROI-1)/2, sum(N.subjects{:,:}));
f = nan(N.TR-1, sum(N.subjects{:,:}));
P = nan(N.TR-1, sum(N.subjects{:,:}));
Xhat = nan(N.ROI*(N.ROI-1)/2, N.TR, sum(N.subjects{:,:}), numel(labels.methods));
msqe = nan(N.TR, N.conditions, size(Phi,4), numel(N.modes)+1);

% Run subject-level DMD
for s = 1:sum(N.subjects{:,:})
    X = FNC{s}(:, 1:N.TR-1);
    Y = FNC{s}(:, 2:N.TR);

    % Run DMD and reconstruction
    for m = 1:numel(labels.methods)
        [Phi(:,:,s,m), mu(:,s), lambda(:,s), diagS(:,s), x0(:,s)] = DMD(X, Y, 'dt',2, 'exact',logical(strcmpi(labels.methods(m),"Exact")), 'r',N.TR-1);
        [Xhat(:,:,s,m), ~] = DMD_recon(Phi(:,:,s,m), lambda(:,m), x0(:,m), N.TR);    % compute reconstruction for each method & number of modes
        msqe(:,s,m) = rmse(Xhat(:,:,s,m), FNC{s});                                   % Compute MSE per sample (TR)
    end

    % Compare exact vs. standard DMD
    D = Phi(:,:,s,1) - Phi(:,:,s,2);
    d = nnz(abs(D) >= eps);
    if d
        dstnc(:,:,s) = abs(Phi(:,:,s,1) - Phi(:,:,s,2));
    end

    % compute true subject-level static FNCs
    sFNC(:,s) = mean(FNC{s}, 2, "omitmissing");

    % Compute group spectra and mode power
    [f(:,s), P(:,s)] = DMD_spectrum(Phi(:,:,s, strcmpi(labels.methods, "exact")), mu(:,s), 'plotit',0);     % power
end

% Check if exact and standard SVD produce same outputs
d = squeeze(sum(dstnc, [1 2]));
if nnz(d)
    warning("DMD methods do not concur!");
end
clear i g s d m n


%% Save results & figure(s)

% % Save figures
% savefig(F, fullfile(pth{5}, fileName), 'compact');
% for c = 1:N.fig-1
%     saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'svg');
%     saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'jpeg');
% end
% clear c F a ax axes ts

% Save files
N.fig = N.fig - 1;
save(fullfile(pth{5}, fileName));
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
fileName = "groups_demeaned";

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
Phi = nan(N.ROI*(N.ROI-1)/2, N.TR-1, N.conditions, numel(labels.methods));
mu = nan(N.TR-1, N.conditions);
lambda = nan(N.TR-1, N.conditions);
diagS = nan(N.ROI*(N.ROI-1)/2, N.conditions);
x0 = nan(N.ROI*(N.ROI-1)/2, N.conditions);
dstnc = cell(N.conditions,1);

% Run group-level DMD
for g = 1:N.conditions
    % Generate group-level X, Y matrices
    X = FNC(analysis_data{:,'Diagnosis'} == labels.diagnosis(g))';
    Y = FNC(analysis_data{:,'Diagnosis'} == labels.diagnosis(g))';
    for s = 1:numel(X)
        X{s} = X{s}(:, 1:N.TR-1);
        Y{s} = Y{s}(:, 2:N.TR);
    end
    X = cell2mat(X);
    Y = cell2mat(Y);

    % Run DMD
    for s = 1:numel(labels.methods)
        [Phi(:,:,g,s), mu(:,g), lambda(:,g), diagS(:,g), x0(:,g)] = DMD(X, Y, 'dt',2, 'exact',logical(strcmpi(labels.methods(s),"Exact")), 'r',N.TR-1);  % standard
    end

    % Compare exact vs. standard DMD
    D = Phi(:,:,g,1) - Phi(:,:,g,2);
    d = nnz(abs(D) >= eps);
    if d
        dstnc{g} = abs(Phi(:,:,g,1) - Phi(:,:,g,2));
    else
        dstnc{g} = [];
    end
end

% compute true group-level static FNCs
sFNC = nan(N.ROI*(N.ROI-1)/2, N.conditions);
for g = 1:N.conditions
    m = cell2mat(FNC(analysis_data{:,"Diagnosis"} == labels.diagnosis(g))');
    sFNC(:,g) = mean(m,2, "omitmissing");
end

% Check if exact and standard SVD produce same outputs
i = cellfun(@isempty, dstnc);
if nnz(i) == length(dstnc)
    Phi = squeeze(Phi(:,:,:,2));    % keep exact DMD
    dstnc = zeros(N.ROI*(N.ROI-1)/2, N.TR-1, N.conditions);
    Xhat = nan(N.ROI*(N.ROI-1)/2, N.TR, N.conditions);

    % Visualize true vs. reconstructed sFNCs
    for g = 1:N.conditions
        F(N.fig) = figure; N.fig = N.fig + 1;
        F(N.fig-1).OuterPosition = [1 1 1920 1055];

        % Plot true group sFNC
        subplot(2,2,1);
        display_FNC(icatb_vec2mat(sFNC(:,g))); hold on
        title(strjoin(["True sFNC of", labels.diagnosis(g)]));

        % Display group dFNC reconstruction
        [Xhat(:,:,g), ~] = DMD_recon(Phi(:,:,g), lambda(:,g), x0(:,g), N.TR-1);
        subplot(2,2,2); title(strjoin(["Reconstructed sFNC of", labels.diagnosis(g)]));
        display_FNC(icatb_vec2mat(mean(Xhat(:,:,g),2))); hold on

        % Plot difference between true, reconstructed sFNC
        d = sFNC(:,g) - mean(Xhat(:,:,g),2);
        subplot(2,2,4); title("True sFNC - Reconstructed sFNC", labels.diagnosis(g));
        display_FNC(icatb_vec2mat(d)); hold on
    end
else
    warning("DMD methods do not concur!");

    % visualize results
    F(N.fig) = figure; N.fig = N.fig + 1;
    F(N.fig-1).OuterPosition = [1 1 1920 1055];
    dstnc = zeros(N.ROI*(N.ROI-1)/2, N.TR-1, N.conditions);
    Xhat = nan(N.ROI*(N.ROI-1)/2, N.TR-1, N.conditions, 2);
    for g = 1:N.conditions

        % Plot true group sFNC
        subplot(2, size(Phi,4)+2, (g-1)*(size(Phi,4)+2)+1);
        display_FNC(icatb_vec2mat(sFNC(:,g)), [0.25 1.5]); hold on
        title(strjoin(["True sFNC of", labels.diagnosis(g)]));

        dstnc(:,:,g) = abs(Phi(:,:,g,1) - Phi(:,:,g,2));

        % reconstruct from DMD
        for s = 1:size(Phi,4)
            [Xhat(:,:,g,s), ~] = DMD_recon(Phi(:,:,g,s), lambda(:,g), x0(:,g), N.TR-1);
            subplot(2, size(Phi,4)+2, (g-1)*(size(Phi,4)+2)+1+s);
            display_FNC(real(icatb_vec2mat(mean(Xhat(:,:,g,s),2))), [0.25 1.5]); hold on;
            title(strjoin(["Reconstructed", labels.diagnosis(g), "sFNC"]), strjoin(["(", lower(labels.methods(s)), ")"], ''));
        end

        % plot difference betwen standard, exact DMD
        subplot(2, size(Phi,4)+2, g*(size(Phi,4)+2));
        display_FNC(real(icatb_vec2mat(mean(Xhat(:,:,g,1)-Xhat(:,:,g,2),2))), [0.25 1.5]);
        title(strjoin(["Reconstructed", labels.diagnosis(g), "sFNC"]), "(standard-exact)"); hold on;
    end
end
clear i g s d m


%% Compute spectra for each group

% preallocate arrays
f = nan(N.TR-1, N.conditions);
P = nan(N.TR-1, N.conditions);

% Compute group spectra and mode power
for g = 1:N.conditions
    % Compute and plot spectra
    [f(:,g), P(:,g), F(N.fig+2*(g-1))] = DMD_spectrum(Phi(:,:,g,1), mu(:,g), 'plotit',1);  % power
    F(N.fig+2*(g-1)).OuterPosition = [1 1 1055 1055]; hold on;   % increase figure size
    title(strjoin(["Power Spectrum for", labels.diagnosis(g)]));
    xlim([min(f(:,g)) max(f(:,g))]); ylim([0 max(P(:,g))]);

    % Plot cumulative power of the modes
    [f_sort, i] = sort(f(:,g),1);
    P_sort = P(i,g);
    F(N.fig+(2*g-1)) = figure; F(N.fig+(2*g-1)).OuterPosition = [1 1 1055 1055];
    plot(f_sort, cumsum(P_sort)./max(cumsum(P_sort),[],'all')); hold on;
    plot([min(f_sort) max(f_sort)], [0.9 0.9], '-r');
    xlabel("frequency (Hz)"); ylabel("% Cumulative Power");
    xlim([min(f_sort) max(f_sort)]);
    title(strjoin(["Cumulative Power for", labels.diagnosis(g), "Group"]));
    legend("Cumulative Power", "90%2 of Power", 'Location','southeast');
end
N.fig = N.fig + 2*g;
clear g X Y d D i P_sort f_sort lambda_sort


%% Check reconstruction error

% get N.modes most powerful modes
i = true(numel(N.modes)+1, N.TR-1);
[~, ind] = sort(P);
for m = 1:numel(N.modes)
    i(m+1, ind(1:N.modes(m),:)) = false;
    i(m+1,:) = ~i(m+1,:);
end
clear ind m

% preallocate for MSQE calculation
mFNC = nan(N.ROI*(N.ROI-1)/2, N.TR, N.conditions);
msqe = nan(N.TR, N.conditions, size(Phi,4), numel(N.modes)+1);
Xhat = nan(N.ROI*(N.ROI-1)/2, N.TR, N.conditions, 2, numel(N.modes)+1);

% Compute and visualize MSQE, reconstructions
for g = 1:N.conditions

    % Compile dFNC for each group
    m = cell2mat(FNC(analysis_data{:,"Diagnosis"} == labels.diagnosis(g))');
    m = reshape(m, N.ROI*(N.ROI-1)/2, N.TR, N.subjects{:,labels.diagnosis(g)});
    mFNC(:,:,g) = mean(m,3);

    for s = 1:numel(labels.methods)     % test both standard and exact DMD

        % Open figure (large)
        F(N.fig) = figure; F(N.fig).OuterPosition = [1 1 1055 1055];
        N.fig = N.fig + 1;
        
        % test reconstruction with several numbers of modes
        for m = 1:numel(N.modes)+1
            % compute reconstruction for each method & number of modes
            [Xhat(:,:,g,s,m), ~] = DMD_recon(Phi(:,:,g,s), lambda(:,g), x0(:,g), N.TR, 'keep_modes',i(m,:));
    
            % Compute MSE per sample (TR)
            msqe(:,g,s,m) = rmse(Xhat(:,:,g,s,m), mFNC(:,:,g));

            % Visualize estimated sFNC
            subplot(2, numel(N.modes)+1, m);
            display_FNC(real(icatb_vec2mat(mean(squeeze(Xhat(:,:,g,s,m)),2))), [0.25 1.5]);
            if m == 1
                title("Reconstructed sFNC (all modes)"); hold on;
            else
                title(strjoin(["Reconstructed sFNC (largest", num2str(nnz(i(m,:))), "modes)"])); hold on;
            end

            % Visualize MSE per sample
            subplot(2, numel(N.modes)+1, m+(numel(N.modes)+1));
            stem(squeeze(msqe(:,g,s,m))); axis tight; hold on
            xlabel('samples'); ylabel('MSE');
            if m == 1
                title("Reconstructed sFNC (all modes)"); hold on;
            else
                title(strjoin(["MSE per sample (largest", num2str(nnz(i(m,:))), "modes)"]));
            end
        end
        sgtitle(F(N.fig-1), strjoin([labels.methods(s), "DMD of", labels.diagnosis(g)]));

        % Open figure (large)
        F(N.fig) = figure; F(N.fig).OuterPosition = [1 1 1055 1055];
        N.fig = N.fig + 1;

        % Visualize group sFNC
        subplot(2,numel(N.modes)+1,1);
        display_FNC(icatb_vec2mat(sFNC(:,g)), [0.25 1.5]); hold on;
        title(strjoin(["sFNC for", labels.diagnosis(g)])); hold on;

        % display difference between reconstructions and true sFNC
        for m = 1:numel(N.modes)+1
            % Display difference between reconstructions and actual sFNC
            subplot(2, numel(N.modes)+1, m+numel(N.modes)+1);
            display_FNC(real(icatb_vec2mat(mean(squeeze(mFNC(:,:,g)) - Xhat(:,:,g,s,m),2))), [0.25 1.5]);
            if m == 1
                title("sFNC - Reconstruction (all modes)"); hold on;
            else
                title(strjoin(["sFNC - Reconstruction (largest", num2str(nnz(i(m,:))), "modes)"])); hold on;
            end
    
            % Display difference between partial and full reconstructions
            if m > 1
                subplot(2, numel(N.modes)+1, m);
                display_FNC(real(icatb_vec2mat(mean(Xhat(:,:,g,s,1) - Xhat(:,:,g,s,m),2))), [0.25 1.5]);
                title("Difference Between Reconstructions", strjoin(["(all modes - largest", num2str(nnz(i(m,:))), "modes)"]));
                hold on;
            end
        end
        sgtitle(F(N.fig-1), strjoin([labels.methods(s), "DMD of", labels.diagnosis(g)]));
    end
end
clear g s m e i mFNC Xhat t


%% Visualize three most powerful modes for each group

% for g = 1:N.conditions
%     for s = 1:numel(labels.methods)
%         for m = 1:numel(N.modes)+1
%             [Phi(:,:,g,s), mu(:,g), lambda(:,g), diagS(:,g), x0(:,g)]
%         end
%     end
% end


%% Plot eigenvalues on unit circle

% test both groups
for g = 1:N.conditions
    % Separate eigenvalues into real, imaginary parts
    i = imag(lambda(:,g));
    r = real(lambda(:,g));
    c(:,1) = abs(lambda(:,g)) > 1;
    c(:,2) = abs(lambda(:,g)) < 1;
    c(:,3) = abs(lambda(:,g)) == 1;
    
    % Plot eigenvalues on unit circle
    F(N.fig) = figure; N.fig = N.fig+1;
    F(N.fig-1).OuterPosition = [1 1 1055 1055];
    pbaspect([1 1 1]); hold on
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
    title(strjoin(labels.diagnosis(g), "Eigenvalues"));
end
clear i r g c theta


%% Plot FN time courses from single module (per mode)

% set mask
r = [6 6 6 7 7 7];
c = [2 3 4 2 3 4];
ind.rc = horzcat(r', c');

% convert masks to linear indices
m = zeros(N.ROI, N.ROI);
m(r,c) = 1;
ind.lin = find(m);

for g = 1:N.conditions                  % test both groups
    [~, i] = sort(f(:,g));              % sort frequencies
    for s = 1:numel(labels.methods)     % test both standard and exact DMD

        % Open figure
        F(N.fig) = figure; N.fig = N.fig+1;
        F(N.fig-1).OuterPosition = [1 1 1920 1055];
        
        % Plot mask
        subplot(2,4,1);
        imagesc(m); colormap bone; colorbar; pbaspect([1 1 1]);
        title("Timecourse Mask");
        xlabel("Neuromark Functional Networks");
        ylabel("Neuromark Functional Networks");
        
        % Plot original FN courses over time
        subplot(2,4,5);
        l = cell2mat(FNC(analysis_data{:,"Diagnosis"} == labels.diagnosis(g))');
        plot(1:N.TR*N.subjects{:,labels.diagnosis(g)}, l(ind.lin(:,1),:)); hold on;
        title("Original FNC Values");
        xlabel("Time Points"); ylabel("Real Amplitude");
        xlim([1 N.TR*N.subjects{:,labels.diagnosis(g)}]);
        legend(num2str(ind.rc));
        
        % Plot FN courses over time as a function of number of modes
        ii = [2 3 4 6 7 8];
        for k = 1:length(c)
            [Xhat, ~] = DMD_recon(Phi(:,i(k+1),g,s), lambda(i(k+1),g), x0(:,g), N.TR*N.subjects{:,labels.diagnosis(g)});    % five most powerful modes
            subplot(2,4,ii(k));
            plot(1:N.TR*N.subjects{:,labels.diagnosis(g)}, Xhat(ind.lin,:)); hold on;
            title("Reconstructed FNC Values", strjoin(["f =", num2str(f(i(k+1))), "Hz"]));
            xlabel("Time Points"); ylabel("Real Amplitude");
            xlim([1 N.TR*10]);
            legend(num2str(ind.rc));
        end

        % title for grid
        sgtitle(strjoin([labels.diagnosis(g), ", ", labels.methods(s), " DMD"], ''));
    end
end
clear i ii ind r c k m n e Xhat l


%% Visualize modes

% test both groups
for g = 1:N.conditions

    % sort frequencies, power
    [f_sort, i] = sort(f(:,g));
    P_sort = P(i,g);

    % test both standard and exact DMD
    for s = 1:numel(labels.methods)

        % remove duplicate (negative) modes
        Phi_sort = Phi(:,i,g,s);
        
        % get amplitudes as function of frequency
        l.r = max(abs(real(Phi_sort)));
        l.i = max(abs(imag(Phi_sort)));
        l.t = max(abs(Phi_sort));

        % visualize dominant harmonic modes
        for j = 1:nnz(cumsum(P_sort)./max(cumsum(P_sort)) < 0.1)
            Phi_mat = icatb_vec2mat(squeeze(Phi_sort(:,j)));
            F(N.fig) = figure; F(N.fig).OuterPosition = [1 1 1920 1055]; N.fig = N.fig + 1;
            subplot(1,2,1); display_FNC(real(Phi_mat), [0.25 1.5]);
            title("Mode (Real Part)"); hold on;
            subplot(1,2,2); display_FNC(imag(Phi_mat), [0.25 1.5], [-max(l.i) max(l.i)]);
            title("Mode (Imaginary Part)"); hold on;
            sgtitle(strjoin([labels.diagnosis(g), ", ", labels.methods(s) " DMD, f = " , num2str(f_sort(j))], ""));
        end
    end

    % visualize amplitudes as function of frequency
    F(N.fig) = figure; F(N.fig).OuterPosition = [1 1 1100 1055]; N.fig = N.fig+1;
    plot(f_sort, l.r, 'r'); hold on
    plot(f_sort, l.i, 'b');
    plot(f_sort, l.t, 'k');
    title('Absolute Amplitudes by Frequency,', labels.diagnosis(g));
    xlabel('Frequency (Hz)'); ylabel('Amplitude');
    legend('Real', 'Imaginary', 'Total');
    xlim([0 max(f,[],"all")]);
end
clear i j Phi_mat phase_mat Phi_sort phi_sort f_sort P_sort l s g


%% Save results & figure(s)

% Save figures
savefig(F, fullfile(pth{5}, fileName), 'compact');
for c = 1:N.fig-1
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'svg');
    saveas(F(c), fullfile(pth{5}, "Images", strjoin([fileName, num2str(c)], '-')), 'jpeg');
end
clear c F a ax axes ts

% Save files
N.fig = N.fig - 1;
save(fullfile(pth{5}, fileName));
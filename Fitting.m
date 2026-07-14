%% D2D
% Historical local fitting script for RON D2D models.
% Current HPC runs use the Setup_RON_*_compile.m and Run_RON_*_LHS.m scripts.
% Configure local paths with environment variables when needed:
%   D2D_FRAMEWORK_DIR, RON_D2D1_DIR, RON_D2D2_DIR, RON_REFERENCE_DIR.
clear all
paths = ron_local_paths();
addpath(genpath(paths.d2d_framework))

% Initialize the Data2Dynamics (D2D) framework
arInit;

cd(paths.d2d1)

% Number of mutant datasets to load per rate parameter category
lenk1 = 62;    % Number of k1 mutants
lenk2 = 274;   % Number of k2 mutants
lenk2 = 320;   % Number of k2 mutants % control
lenk3 = 69;    % Number of k3 mutants

% Choose model type: 1 for one-step model, else two-step model
model_type = 3;


% Load appropriate model based on model_type
if model_type == 1
    arLoadModel('one_step');
    data_label = 'RON_HEK293';
elseif model_type==2;
    arLoadModel('two_step');
    data_label = 'RON_HEK293';
elseif model_type==3;%control fitting
    arLoadModel('two_step');
    data_label = 'RON_control';
end

% Define experiment/dataset label

% Load WT (wild type) dataset
arLoadData(data_label);

% Load mutant datasets for k1 variants
% for i = 1:lenk1
%     arLoadData([data_label '_mut_k1_' num2str(i)]);
% end

% Load mutant datasets for k2 variants
for i = 1:lenk2
    arLoadData([data_label '_mut_k2_' num2str(i)]);
end

% Load mutant datasets for k3 variants
% for i = 1:lenk3
%     arLoadData([data_label '_mut_k3_' num2str(i)]);
% end

% Compile the model with all loaded data
arCompileAll;

% Define parameter search boundaries
limkon  = [-2 -3 1];    % log10 bounds for kon rates (association)
limkoff = [-2];         % initial value for koff (dissociation)
limdeg  = [-5 -6 -4];   % log10 bounds for degradation rates
limkspli = [-1 -3 1];   % log10 bounds for splicing rates

% Set kinetic parameters with log10-transformed bounds and initial values
arSetPars('k1', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k2', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k2a', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k2b', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k3', limkon(1), 1, 1, limkon(2), limkon(3));

arSetPars('kret', limkoff, 0, 1, -4, 1);  % return rate, fixed

% koff variants
arSetPars('k4',  limkoff, 0, 1, -5, 5);
arSetPars('k5a', limkoff, 0, 1, -5, 5);
arSetPars('k5b', limkoff, 0, 1, -5, 5);
arSetPars('k5',  limkoff, 0, 1, -5, 5);
arSetPars('k6',  limkoff, 0, 1, -5, 5);

% Splicing rate
arSetPars('kspli', limkspli(1), 1, 1, limkspli(2), limkspli(3));

% Degradation rates
arSetPars('kincl',  limdeg(1), 1, 1, limdeg(2), limdeg(3));  % Inclusion
arSetPars('kskip',  limdeg(1), 1, 1, limdeg(2), limdeg(3));  % Skipping
arSetPars('kdr1',   limdeg(1), 1, 1, limdeg(2), limdeg(3));  % Degradation 1
arSetPars('kdr2',   limdeg(1), 1, 1, limdeg(2), limdeg(3));  % Degradation 2

% Scaling parameter (e.g., signal strength), bounded
arSetPars('s', -3, 0, 1, -5, -2);


% arLoad('Results/20251003T165249_two_step_ready');

% Disable fitting of error parameters
ar.config.fiterrors = 0;

% Perform model fitting using Latin Hypercube Sampling
fit = 1000; 
arFitLHS(fit)

% Save results to a .mat file
model_type = 3;

if model_type == 1
    data_label = 'RON_HEK293';
elseif model_type==2;
    data_label = 'RON_HEK293';
elseif model_type==3;%control fitting
    data_label = 'RON_control';
end

% Save results to a .mat file
filename = [data_label, '_', num2str(model_type), '_', num2str(fit),'_SE',  '.mat'];
save(filename);

ar.config.saveAllPars = 1;   
saveDir = [ar.model.name, '_', num2str(fit), '_SE'];
arSave(saveDir);


%% Using Pananjots parameter
% === SETUP ===
clear all
paths = ron_local_paths();
addpath(genpath(paths.d2d_framework))
arInit;
% ar.config.checksum = false;

cd(paths.d2d2)

% === Load compiled model and data ===
model_type = 1;  % 1 = one-step; 2 = two-step; 3 = control
if model_type == 1
    arLoadModel('one_step');
    data_label = 'RON_HEK293';
elseif model_type == 2
    arLoadModel('two_step');
    data_label = 'RON_HEK293';
elseif model_type == 3
    arLoadModel('two_step');
    data_label = 'RON_control';
end

arLoadData(data_label);
lenk1 = 62; lenk2 = 274; lenk3 = 69;
for i = 1:lenk1
    arLoadData([data_label '_mut_k1_' num2str(i)]);
end
for i = 1:lenk2
    arLoadData([data_label '_mut_k2_' num2str(i)]);
end
for i = 1:lenk3
    arLoadData([data_label '_mut_k3_' num2str(i)]);
end

ar.config.nParallelCompile = 8; 
arCompileAll;

cd(paths.reference_results)

% parameters corresponding to the best fit result
S = load('RON_HEK293_1_1000.mat', 'ar'); 
old_ar = S.ar;     
[best_chi2, best_idx] = min(old_ar.chi2s); 
best_pars = old_ar.ps(best_idx, :)'; 

% Copy old boundaries into the current model
for i = 1:length(ar.pLabel)
    lbl = ar.pLabel{i};
    idx_old = find(strcmp(old_ar.pLabel, lbl));
    if ~isempty(idx_old)
        ar.lb(i) = old_ar.lb(idx_old);
        ar.ub(i) = old_ar.ub(idx_old);
    end
end

% Set it to the current model
for i = 1:length(ar.p)
    ar.p(i) = best_pars(i);
end

cd(paths.d2d2)
ar.config.fiterrors = 0;
% === Single-point fitting from best guess ===
arFit;
ar.config.saveAllPars = 1;
arSave('one_step_refined_from_LHS_best')
% 
% % Save results to a .mat file
% filename = [data_label, '_', num2str(model_type), '_PanajotV1','.mat'];
% save(filename);

% Perform model fitting using Latin Hypercube Sampling
fit = 1000; 
arFitLHS(fit);
% === Save results ===
ar.config.saveAllPars = 1;
arSave('one_step_refined_from_LHS_best_1000')
% Save results to a .mat file
filename = [data_label, '_', num2str(model_type), '_', num2str(fit),'_PanajotV1', '.mat'];
save(filename);


%% Check
%%%%%%%%%%%  lower bound & upper bound  %%%%%%%%%%%
% parameter labels
labels = old_ar.pLabel(:);
% best parameters from Panajot
pars = best_pars(:);
% lower & upper bounds
lb = old_ar.lb(:);
ub = old_ar.ub(:);
% Assemble table
T = table(labels, pars, lb, ub);
disp(T);

idx_low  = pars < lb;
idx_high = pars > ub;
fprintf('\n=== Parameters BELOW lower bound ===\n');
disp(T(idx_low, :));
fprintf('\n=== Parameters ABOVE upper bound ===\n');
disp(T(idx_high, :));



%%%%%%%%%%%  arSimu check  %%%%%%%%%%%
%%% arCompileAll %%%
clear all
paths = ron_local_paths();
addpath(genpath(paths.d2d_framework))
arInit;
%  Load model and data (same as Panajot used) 
cd(paths.reference_results)
arLoadModel('one_step');
data_label = 'RON_HEK293';
arLoadData(data_label);
lenk1 = 62; lenk2 = 274; lenk3 = 69;
for i = 1:lenk1
    arLoadData([data_label '_mut_k1_' num2str(i)]);
end
for i = 1:lenk2
    arLoadData([data_label '_mut_k2_' num2str(i)]);
end
for i = 1:lenk3
    arLoadData([data_label '_mut_k3_' num2str(i)]);
end
%  RECOMPILE model (this regenerates the missing MEX file)
arCompileAll;

%%% arSimu %%% 
% clear all
% addpath(genpath(paths.d2d_framework))
% arInit;

S = load(fullfile(paths.reference_results, 'RON_Hek293_1_1000.mat'));
old_ar = S.ar;
[~, best_idx] = min(old_ar.chi2s);
best_pars = old_ar.ps(best_idx, :)';
% apply Panajot parameters to current model
ar.p = best_pars;
% Now simulate 
arSimu;

% Collect all predicted isoforms
pred = vertcat(ar.model.data.yExpSimu);
% Make PSI–Efficiency plot
psi = pred(:,1)./(pred(:,1)+pred(:,2));
eff = 1 - sum(pred(:,3:5),2);
figure; plot(psi*100, eff*100, '.');
axis square
xlabel('PSI (%)'); ylabel('Efficiency (%)');
title('PSI–Efficiency (Panajot original one-step)');

% Make PSI–Efficiency plot -- Just k2 mutants
% Define index for k2 mutants
idx_k2 = (lenk1 + 2) : (lenk1 + lenk2 + 1);
pred_k2 = pred(idx_k2, :);
% Make PSI–Efficiency plot
psi = pred_k2(:,1) ./ (pred_k2(:,1)+pred_k2(:,2));
eff = 1 - sum(pred_k2(:,3:5),2);
figure; 
plot(psi*100, eff*100, '.');
axis square
xlabel('PSI (%)'); ylabel('Efficiency (%)');
title('PSI–Efficiency (k2 mutants only)');

%%%%%%%%%%% Multi-start arFit using Panajot best parameters %%%%%%%%%%%
nFits = 2;
chi2_hist = zeros(nFits,1);
psi_runs = cell(nFits,1);
eff_runs = cell(nFits,1);

for k = 1:nFits
    ar.p = best_pars;          % restart from best point every run
    arFit;                     % single-point refinement
    chi2_hist(k) = arGetMerit('chi2fit');
    
    ySimAll = vertcat(ar.model(1).data.yExpSimu);   % predicted means
    pred_k2 = ySimAll((lenk1 + 2):(lenk1 + lenk2 + 1), :);
    psi_runs{k} = pred_k2(:,1) ./ (pred_k2(:,1)+pred_k2(:,2));
    eff_runs{k} = 1 - sum(pred_k2(:,3:5),2);
end

figure;
hold on;
for k = 1:nFits
    plot(psi_runs{k}*100, eff_runs{k}*100, '.');
end
axis square
xlabel('PSI (%)');
ylabel('Efficiency (%)');
title('PSI–Efficiency (k2 mutants only)');
hold off;

save('Repeated_arFit_k2_2.mat', ...
     'best_pars', ...
     'chi2_hist', ...
     'psi_runs', ...
     'eff_runs', ...
     'lenk1','lenk2','lenk3');


%%%% V2 %%%% 
Nrepeat = 10;                 
chi2_list = zeros(Nrepeat,1);
psi_all = cell(Nrepeat,1);
eff_all = cell(Nrepeat,1);
pars_all = zeros(Nrepeat, length(ar.p));

ar.p = best_pars; 
for r = 1:Nrepeat
    % --- Run local fit ---
    arFit;
    chi2_list(r) = ar.chi2;
    pars_all(r,:) = ar.p';

    % --- Simulate with fitted parameters ---
    arSimu;

    % --- Collect predicted isoforms ---
    pred = vertcat(ar.model.data.yExpSimu);

    % --- Extract PSI / Efficiency (all mutants or subset) ---
    idx_k2 = (lenk1 + 2) : (lenk1 + lenk2 + 1);
    psi = pred(idx_k2,1) ./ (pred(idx_k2,1) + pred(idx_k2,2));
    eff = 1 - sum(pred(idx_k2,3:5), 2);

    psi_all{r} = psi;
    eff_all{r} = eff;
end

% Plot PSI–Efficiency curves overlaying all runs 
figure; hold on;
for r = 1:Nrepeat
    plot(psi_all{r}*100, eff_all{r}*100, '.', 'MarkerSize', 8);
end
xlabel('PSI (%)'); ylabel('Efficiency (%)');
title(sprintf('PSI–Efficiency overlay for %d repeated arFit runs', Nrepeat));
axis square; grid on;

%  Inspect chi2 across runs 
figure; plot(chi2_list,'o-','LineWidth',2);
xlabel('Run number'); ylabel('chi^2');
title('chi^2 across repeated arFit runs');

% Save all results 
save('Repeated_arFit_results_2.mat', ...
     'chi2_list','pars_all','psi_all','eff_all');

%% Get & save Panajots Parameter
% Add all subfolders from 'data2d-master' to MATLAB path
% addpath(genpath(paths.d2d_framework))

arInit;
load('RON_Hek293_1_1000.mat', 'ar');   

T = table( ...
    ar.p(:), ...
    ar.lb(:), ...
    ar.ub(:), ...
    ar.qFit(:), ...
    ar.pLabel(:), ...
    'VariableNames', {'p','lb','ub','qFit','pLabel'});

writetable(T,'pars_panajot_1.txt');



%% Start Pananjots parameter
arInit;
% CD 
paths = ron_local_paths();
cd(paths.d2d1)

lenk1 = 62;    
lenk2 = 274;   
lenk3 = 69;    
model_type = 3; 

if model_type == 1
    arLoadModel('one_step');
    data_label = 'RON_HEK293';
elseif model_type==2
    arLoadModel('two_step');
    data_label = 'RON_HEK293';
elseif model_type==3
    arLoadModel('two_step');
    data_label = 'RON_control';
end

arLoadData(data_label);
for i = 1:lenk1; arLoadData([data_label '_mut_k1_' num2str(i)]); end
for i = 1:lenk2; arLoadData([data_label '_mut_k2_' num2str(i)]); end
for i = 1:lenk3; arLoadData([data_label '_mut_k3_' num2str(i)]); end

arCompileAll;

cd(paths.reference_results)


% parameters corresponding to the best fit result
S = load('RON_HEK293_1_1000.mat', 'ar'); 
old_ar = S.ar;     
[best_chi2, best_idx] = min(old_ar.chi2s); 
best_pars = old_ar.ps(best_idx, :)'; 

% --- map EVERYTHING by label (p, lb, ub, qFit, qLog10) ---
for i = 1:length(ar.pLabel)
    lbl = ar.pLabel{i};
    idx_old = find(strcmp(old_ar.pLabel, lbl), 1);
    if ~isempty(idx_old)
        ar.p(i)    = best_pars(idx_old);   % <-- KEY: by label
        ar.lb(i)   = old_ar.lb(idx_old);
        ar.ub(i)   = old_ar.ub(idx_old);
        ar.qFit(i) = old_ar.qFit(idx_old);

        % only if both sides have qLog10
        if isfield(old_ar,'qLog10') && isfield(ar,'qLog10') && ~isempty(old_ar.qLog10)
            ar.qLog10(i) = old_ar.qLog10(idx_old);
        end
    else
        warning('No match for parameter label: %s', lbl);
    end
end


try
    [ar, ~] = arCalcMerit(ar, false, ar.p(ar.qFit==1));  % no sensis
    fprintf('chi2=%.6g, anyNaNres=%d, anyInfres=%d\n', ar.chi2, any(isnan(ar.res)), any(isinf(ar.res)));
catch ME
    disp(ME.message);
end



% fprintf(' Enable all parameters (ar.qFit) for fitting....\n');
% ar.qFit = ones(size(ar.qFit));

cd(paths.d2d2)
ar.config.fiterrors = 0;
n_runs = 500;
perturb_range = 0.5;

% [ar, merit] = arCalcMerit(ar, false, ar.p(ar.qFit == 1));

arFit;

best_chi2 = ar.chi2;
best_ar   = ar;

% chi2 history
chi2_history = zeros(n_runs, 1);     
best_chi2_over_time = zeros(n_runs, 1); 

for i = 1:n_runs
    fprintf('\n=== Multistart run %d / %d ===\n', i, n_runs);
    tStart = tic;

    ar.p = best_ar.p; 
    
    % small perturbations  --- only to the fitable parameters (ar.qFit == 1).
    rand_perturbation = (rand(size(ar.p)) * 2 - 1) * perturb_range; 
    ar.p(ar.qFit == 1) = ar.p(ar.qFit == 1) + rand_perturbation(ar.qFit == 1);
   
    % do not exceed the boundary.
    ar.p = min(ar.ub, max(ar.lb, ar.p));

    arFit; 

    elapsed = toc(tStart);
    
    % Compare and save the results
    current_chi2 = ar.chi2;
    chi2_history(i) = current_chi2;

    if current_chi2 < best_chi2
        fprintf('  A better solution has been found!！Chi2: %f -> %f\n', best_chi2, current_chi2);
        best_chi2 = current_chi2;
        best_ar = ar; 
    else
        fprintf('  ❌ The result did not improve.Chi2: %f\n', current_chi2);
    end
    best_chi2_over_time(i) = best_chi2;
end

% save
ar = best_ar;

%% Plot
% Chi2 History Plot
figure;
hold on;
plot_range = 1:n_runs;
% chi2_history
h1 = plot(plot_range, chi2_history, 'o', 'MarkerSize', 3, 'MarkerEdgeColor', [0.5 0.5 0.5], 'DisplayName', 'Final Chi^2 of Each Run');
% best_chi2_over_time
h2 = plot(plot_range, best_chi2_over_time, 'r-', 'LineWidth', 2, 'DisplayName', 'Best Chi^2 Found So Far');
% Mark the final optimal Chi2 
[final_best, best_run_idx] = min(chi2_history);
plot(best_run_idx, final_best, 'b*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Global Best Run');
grid on;
box on;
xlabel('Run Number');
ylabel('Chi^2 Value');
title(sprintf('Optimization History: Chi^2 over %d Multi-start Runs', n_runs));
legend('show');
hold off;


% PSI–Efficiency Plot
arCalcRes;
pred = vertcat(ar.model.data.yExpSimu);
% PSI = Inclusion / (Inclusion + Skipping)
psi = pred(:,1) ./ (pred(:,1) + pred(:,2));
eff = 1 - sum(pred(:,3:5), 2);
figure; 
plot(psi*100, eff*100, 'o', 'MarkerSize', 5, 'DisplayName', 'All Mutants'); 
axis square;
grid on;
xlabel('PSI (%)'); 
ylabel('Efficiency (%)');
title('PSI–Efficiency Plot (Refitted Model)');

% ---  k2  ---
idx_k2_start = lenk1 + 2; 
idx_k2_end = lenk1 + lenk2 + 1;
idx_k2 = idx_k2_start : idx_k2_end; 

if idx_k2_end <= size(pred, 1)
    pred_k2 = pred(idx_k2, :);
    
    psi_k2 = pred_k2(:,1) ./ (pred_k2(:,1) + pred_k2(:,2));
    eff_k2 = 1 - sum(pred_k2(:,3:5), 2);
    
    figure; 
    plot(psi_k2*100, eff_k2*100, 'r.', 'MarkerSize', 10);
    axis square;
    grid on;
    xlabel('PSI (%)'); 
    ylabel('Efficiency (%)');
    title('PSI–Efficiency (k2 mutants only)');
end

% Model vs Data Plot
arCalcRes;
dr = ar.model.data;
data_all = vertcat(dr.yExp);
pred_all = vertcat(dr.yExpSimu);

[N, M] = size(data_all);
isoNames = {'Inclusion', 'Skipping', 'FULL IR', '1IR', '2IR'}; 

figure;
hold on;
% colors and markers
colorOrder = get(gca, 'ColorOrder');
h_plots = []; 
for m = 1:M 
    h_plots(m) = plot(data_all(:, m) * 100, pred_all(:, m) * 100, 'o', ...
        'Color', colorOrder(mod(m-1, size(colorOrder, 1)) + 1, :), ...
        'MarkerFaceColor', colorOrder(mod(m-1, size(colorOrder, 1)) + 1, :), ...
        'MarkerSize', 4, ...
        'DisplayName', isoNames{m});
end
% Bisecting Line 
plot([0 100], [0 100], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Perfect Fit'); 
% Model error area using previously fitted curve parameters p.
% old p without refitting p may introduce errors.
if exist('p', 'var') 
    yw = linspace(0, 1);
    vv = polyval(p, yw);
    fill([yw, fliplr(yw)] * 100, [yw - vv, fliplr(yw + vv)] * 100, ...
        [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'LineWidth', 0.1, 'HandleVisibility', 'off');
end
xlabel('Data (%)');
ylabel('Model (%)');
xlim([0 100]);
ylim([0 100]);
title('Model Prediction vs Experimental Data');
legend('show', 'Location', 'SouthEast');
axis square;
hold off;

function paths = ron_local_paths()
%RON_LOCAL_PATHS Resolve local D2D folders without hard-coded user paths.
% Override these defaults with environment variables before running sections:
%   D2D_FRAMEWORK_DIR, RON_D2D1_DIR, RON_D2D2_DIR, RON_REFERENCE_DIR.

home_dir = getenv('HOME');
paths.d2d_framework = getenv_or_default( ...
    'D2D_FRAMEWORK_DIR', fullfile(home_dir, 'HPC', 'd2d-master', 'arFramework3'));
paths.d2d1 = getenv_or_default( ...
    'RON_D2D1_DIR', fullfile(home_dir, 'HPC', 'D2D_1'));
paths.d2d2 = getenv_or_default( ...
    'RON_D2D2_DIR', fullfile(home_dir, 'HPC', 'D2D_2'));
paths.reference_results = getenv_or_default( ...
    'RON_REFERENCE_DIR', fullfile(home_dir, 'HPC', 'reference_results'));
end

function value = getenv_or_default(name, default_value)
value = getenv(name);
if isempty(value)
    value = default_value;
end
end










%% Run_RON_onestep_LHS.m
% Run one-step RON model fitting as a SLURM array job.
% Requires RESULT_DIR and optionally SLURM_ARRAY_TASK_ID from the scheduler.
clear all; close all;

%% --- paths ---
home_dir = getenv('HOME');
proj_root = fullfile(home_dir, 'HPC');
d2d_dir   = fullfile(proj_root, 'D2D_1');

result_dir = getenv('RESULT_DIR');
if isempty(result_dir)
    error('RESULT_DIR not set (did you run via sbatch?)');
end

addpath(genpath(fullfile(proj_root, 'd2d-master', 'arFramework3')));
cd(d2d_dir);

%% --- init D2D ---
arInit;
ar.config.useSensis    = 0;
ar.config.useSensisErr = 0;
ar.config.fiterrors    = 0;

%% --- load compiled model ---
load(fullfile(proj_root, 'RON_onestep_compiled_test.mat'), 'ar');

%% --- load old result ONLY for bounds ---
S = load(fullfile(proj_root, 'RON_Hek293_1_1000.mat'), 'ar');
old_ar = S.ar;

for i = 1:length(ar.pLabel)
    idx_old = find(strcmp(old_ar.pLabel, ar.pLabel{i}), 1);
    if ~isempty(idx_old)
        ar.lb(i)   = old_ar.lb(idx_old);
        ar.ub(i)   = old_ar.ub(idx_old);
        ar.qFit(i) = old_ar.qFit(idx_old);
    end
end

%% --- free kret ---
kret_idx = find(strcmp(ar.pLabel, 'kret'));
if ~isempty(kret_idx)
    ar.qFit(kret_idx) = 1;
end

%% --- SLURM array info ---
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(task_id)
    task_id = 1;
end

fit_per_job = 100;
rng(2000 + task_id);

fprintf('=== ONE-STEP LHS batch %d / %d samples ===\n', task_id, fit_per_job);
fprintf('Result dir: %s\n', result_dir);

%% --- run LHS ---
arFitLHS(fit_per_job);

%% --- save ---
outfile = fullfile(result_dir, ...
    sprintf('RON_onestep_LHS_batch_%02d.mat', task_id));

save(outfile, 'ar');

fprintf('Saved results to %s\n', outfile);
fprintf('ONE-STEP LHS batch %d finished.\n', task_id);

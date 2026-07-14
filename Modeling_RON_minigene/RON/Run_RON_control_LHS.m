%% Run_RON_control_LHS.m
% Run two-step control RON model fitting as a SLURM array job.
% Requires RESULT_DIR and optionally SLURM_ARRAY_TASK_ID from the scheduler.
clear; close all;

% === Paths ===
home_dir   = getenv('HOME');
proj_root  = fullfile(home_dir, 'HPC');
d2d_dir    = fullfile(proj_root, 'D2D_1');
result_dir = getenv('RESULT_DIR');
if isempty(result_dir)
    error('RESULT_DIR not set (did you run via sbatch?)');
end

addpath(genpath(fullfile(proj_root, 'd2d-master', 'arFramework3')));
cd(d2d_dir);
arInit;

% === Load compiled model ===
load(fullfile(proj_root, 'RON_control_twostep_compiled_SE.mat'), 'ar');

% === Set parameter bounds manually ===
limkon   = [-2 -3 1];
limkoff  = [-2];
limdeg   = [-5 -6 -4];
limkspli = [-1 -3 1];

arSetPars('k1', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k2a', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k2b', limkon(1), 1, 1, limkon(2), limkon(3));
arSetPars('k3', limkon(1), 1, 1, limkon(2), limkon(3));

arSetPars('kret', limkoff, 0, 1, -4, 1);

arSetPars('k4',  limkoff, 0, 1, -5, 5);
arSetPars('k5a', limkoff, 0, 1, -5, 5);
arSetPars('k5b', limkoff, 0, 1, -5, 5);
arSetPars('k6',  limkoff, 0, 1, -5, 5);

arSetPars('kspli', limkspli(1), 1, 1, limkspli(2), limkspli(3));

arSetPars('kincl',  limdeg(1), 1, 1, limdeg(2), limdeg(3));
arSetPars('kskip',  limdeg(1), 1, 1, limdeg(2), limdeg(3));
arSetPars('kdr1',   limdeg(1), 1, 1, limdeg(2), limdeg(3));
arSetPars('kdr2',   limdeg(1), 1, 1, limdeg(2), limdeg(3));

arSetPars('s', -3, 0, 1, -5, -2);

% === Fitting settings ===
ar.config.fiterrors = 0;

% === SLURM task info ===
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(task_id)
    task_id = 1;
end

rng(1000 + task_id);
fit_per_job = 100;

fprintf('=== SLURM Task %d: Fitting %d LHS samples ===\n', task_id, fit_per_job);
fprintf('Result dir: %s\n', result_dir);

arFitLHS(fit_per_job);

% === Save results ===
outfile = fullfile(result_dir, sprintf('RON_control_LHS_batch_%02d.mat', task_id));
save(outfile, 'ar');

fprintf('Saved result to %s\n', outfile);

% Setup_RON_control_compile.m
% Compile the two-step control RON D2D model on the HPC/Linux environment.
% Expected project layout: $HOME/HPC/d2d-master/arFramework3 and $HOME/HPC/D2D_1.
clear all; close all;

home_dir = getenv('HOME');
proj_root = fullfile(home_dir, 'HPC');
d2d_path  = fullfile(proj_root, 'd2d-master', 'arFramework3');
d2d_dir   = fullfile(proj_root, 'D2D_1');

addpath(genpath(d2d_path));
arInit;

cd(d2d_dir);
arLoadModel('two_step');
arLoadData('RON_control');

% Load all mutation data
lenk1 = 115;
lenk2 = 320;
lenk3 = 139;

% for i = 1:lenk1, arLoadData(['RON_control_mut_k1_' num2str(i)]); end

for i = 1:lenk2, arLoadData(['RON_control_mut_k2_' num2str(i)]); end

% for i = 1:lenk3, arLoadData(['RON_control_mut_k3_' num2str(i)]); end

ar.config.nParallelCompile = 8;
arCompileAll;

save(fullfile(proj_root, 'RON_control_twostep_compiled_SE.mat'), 'ar');

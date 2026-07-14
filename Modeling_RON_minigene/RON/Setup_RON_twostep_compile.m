% Setup_RON_twostep_compile.m
% Compile the two-step RON D2D model on the HPC/Linux environment.
% Expected project layout: $HOME/HPC/d2d-master/arFramework3 and $HOME/HPC/D2D_1.
clear all
close all

home_dir = getenv('HOME');

% --- paths inside the HPC project folder ---
proj_root = fullfile(home_dir, 'HPC');     
d2d_path  = fullfile(proj_root, 'd2d-master', 'arFramework3');
d2d_dir   = fullfile(proj_root, 'D2D_1');
 
% --- add D2D ---
addpath(genpath(d2d_path));
arInit;

% --- switch to D2D project directory ---
cd(d2d_dir);

% --- load model & data ---
lenk1 = 62;    
lenk2 = 274;   
lenk3 = 69;    
model_type = 2;   % two-step

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
for i = 1:lenk1, arLoadData([data_label '_mut_k1_' num2str(i)]); end
for i = 1:lenk2, arLoadData([data_label '_mut_k2_' num2str(i)]); end
for i = 1:lenk3, arLoadData([data_label '_mut_k3_' num2str(i)]); end

% --- compile to Linux MEX ---
ar.config.nParallelCompile = 8;
arCompileAll;   

% --- save compiled model object ---
save(fullfile(proj_root, 'RON_twostep_compiled.mat'), 'ar');

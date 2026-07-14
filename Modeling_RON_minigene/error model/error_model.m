%% error_model_Jay.m
% Estimate empirical isoform error models from HEK293 replicate data.
% The main output is error_model_hek_withGLOBAL_allExons_withOffset.mat,
% containing per-isoform and global sigma(mu) functions used by D2D input
% generation.

%% K-transform error model for exon 2
load("datahekrep1nowt.mat");
load("datahekrep2nowt.mat");
load("datahekrep3nowt.mat");

datahekrep1 = datahekrep1nowt; 
datahekrep2 = datahekrep2nowt; 
datahekrep3 = datahekrep3nowt;

% Exon 2
isoCols = 2:6;                      % Inclusion, Skipping, FullIR, FirstIR, SecondIR
[C2, idx2_r1, idx2_r2, idx2_r3] = mintersect( ...
    datahekrep1{3}(:,1), datahekrep2{3}(:,1), datahekrep3{3}(:,1));

R1 = datahekrep1{3}(idx2_r1, isoCols);   % N x 5
R2 = datahekrep2{3}(idx2_r2, isoCols);   % N x 5
R3 = datahekrep3{3}(idx2_r3, isoCols);   % N x 5

X  = cat(3, R1, R2, R3);                 % N x 5 x 3
mu_ex2  = mean(X, 3);                    % N x 5
std_ex2 = std(X, 0, 3);                  % N x 5

isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};




%%% --- 1 --- Fitting by isoform
Kiso = numel(isoCols);
ERR = struct('name',[],'intercept',[],'alpha',[],'model',[],'sigma_fun',[]);

for k = 1:Kiso
    x = mu_ex2(:,k); 
    y = std_ex2(:,k);

    % 1) Linear fit mean vs std, take the intercept
    p_lin = polyfit(x, y, 1);
    b0 = p_lin(2);             % intercept_k

    % 2) K transformation
    K_vals = (1 ./ x) - 1;

    % 3) lsqnonlin fit alpha_k
    errfun = @(alpha) -y + ((1./(1 + K_vals)) - (1./(1 + alpha .* K_vals)) + b0);
    opts = optimoptions('lsqnonlin','Display','off');
    alpha0 = 0;
    alpha_k = lsqnonlin(errfun, alpha0, [], [], opts);

    % 4) σ(μ) function of the isoform
    sigma_fun_k = @(muq) max( ...
        (1./(1 + (1./muq - 1)) - 1./(1 + alpha_k.*(1./muq - 1))) + b0, ...
        0);

    % Save
    ERR(k).name = isoNames{k};
    ERR(k).intercept = b0;
    ERR(k).alpha = alpha_k;
    ERR(k).model = 'std(mu) = f_alpha,b0(mu)'; 
    ERR(k).sigma_fun = sigma_fun_k;
end

%%% --- 1A --- GLOBAL model (all isoforms together)
x_all = mu_ex2(:);
y_all = std_ex2(:);

% 1) Linear fit
p_lin_all = polyfit(x_all, y_all, 1);
b0_all = p_lin_all(2);

% 2) K transformation
K_all = (1 ./ x_all) - 1;

% 3) Fit alpha_all
errfun_all = @(alpha) -y_all + ((1./(1 + K_all)) - (1./(1 + alpha .* K_all)) + b0_all);
opts = optimoptions('lsqnonlin','Display','off');
alpha0 = 0;
alpha_all = lsqnonlin(errfun_all, alpha0, [], [], opts);

% 4) σ(μ) function
sigma_fun_all = @(muq) max( ...
    (1./(1 + (1./muq - 1)) - 1./(1 + alpha_all.*(1./muq - 1))) + b0_all, ...
    0);

% Save into ERR(6)
ERR(Kiso+1).name = 'GLOBAL';
ERR(Kiso+1).intercept = b0_all;
ERR(Kiso+1).alpha = alpha_all;
ERR(Kiso+1).model = 'std(mu) = f_alpha,b0(mu) [global]';
ERR(Kiso+1).sigma_fun = sigma_fun_all;

%%% --- 2 --- Plot all fits
figure; 
tl = tiledlayout(3,2,'TileSpacing','compact'); 
for k = 1:Kiso+1
    nexttile; 
    if k <= Kiso
        x = mu_ex2(:,k); 
        y = std_ex2(:,k);
    else
        x = mu_ex2(:);  
        y = std_ex2(:);
    end
    xfit = linspace(0,1,200);
    yfit = ERR(k).sigma_fun(xfit);
    plot(x*100, y*100, 'o'); hold on;
    plot(xfit*100, yfit*100, 'k','LineWidth',1.3); grid on
    title(sprintf('%s | \\alpha=%.3g, b0=%.3g', ERR(k).name, ERR(k).alpha, ERR(k).intercept));
    xlabel('Mean (%)'); ylabel('Std (%)');
end


%%  Data MEAN 
load("datahekrep1nowt.mat");
load("datahekrep2nowt.mat");
load("datahekrep3nowt.mat");

datahekrep1 = datahekrep1nowt; 
datahekrep2 = datahekrep2nowt; 
datahekrep3 = datahekrep3nowt;

% Align three replicates on exon2 (cell {3}) by exact position intersection
isoCols = 2:6;  % Inclusion, Skipping, FullIR, FirstIR, SecondIR
[C2, idx2_r1, idx2_r2, idx2_r3] = mintersect( ...
    datahekrep1{3}(:,1), datahekrep2{3}(:,1), datahekrep3{3}(:,1));

R1 = datahekrep1{3}(idx2_r1, isoCols);   % N x 5
R2 = datahekrep2{3}(idx2_r2, isoCols);   % N x 5
R3 = datahekrep3{3}(idx2_r3, isoCols);   % N x 5

X  = cat(3, R1, R2, R3);                 % N x 5 x 3
mu_ex2  = mean(X, 3);                    % N x 5
std_ex2 = std(X, 0, 3);                  % N x 5

isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};



%%% --- 1 --- Per-isoform error models (mean-based formula)
Kiso = numel(isoCols);
ERR = struct('name',[],'intercept',[],'alpha',[],'model',[],'sigma_fun',[]);

for k = 1:Kiso
    x = mu_ex2(:,k);     % mean
    y = std_ex2(:,k);    % std

    % 1) Linear fit y = m*x + b, take intercept b0
    p_lin = polyfit(x, y, 1);
    b0 = p_lin(2);

    % 2) Fit alpha by nonlinear least squares using mean-based formula
    %    std(mu) = 1/(1+mu) - 1/(1+alpha*mu) + b0
    errfun = @(alpha) -y + ( 1./(1 + x) - 1./(1 + alpha.*x) + b0 );
    opts = optimoptions('lsqnonlin','Display','off');
    alpha0 = 0;
    alpha_k = lsqnonlin(errfun, alpha0, [], [], opts);

    % 3) σ(μ) function for this isoform (returns std, not variance)
    sigma_fun_k = @(muq) max( 1./(1 + muq) - 1./(1 + alpha_k.*muq) + b0 , 0 );

    % Save
    ERR(k).name = isoNames{k};
    ERR(k).intercept = b0;
    ERR(k).alpha = alpha_k;
    ERR(k).model = 'std(mu) = 1/(1+mu) - 1/(1+alpha*mu) + intercept';
    ERR(k).sigma_fun = sigma_fun_k;
end

%%% --- 1A --- GLOBAL model (all isoforms pooled)
x_all = mu_ex2(:);
y_all = std_ex2(:);

% 1) Global intercept
p_lin_all = polyfit(x_all, y_all, 1);
b0_all = p_lin_all(2);

% 2) Fit global alpha (mean-based)
errfun_all = @(alpha) -y_all + ( 1./(1 + x_all) - 1./(1 + alpha.*x_all) + b0_all );
opts = optimoptions('lsqnonlin','Display','off');
alpha0 = 0;
alpha_all = lsqnonlin(errfun_all, alpha0, [], [], opts);

% 3) Global σ(μ)
sigma_fun_all = @(muq) max( 1./(1 + muq) - 1./(1 + alpha_all.*muq) + b0_all , 0 );

% Store as ERR(6)
ERR(Kiso+1).name = 'GLOBAL';
ERR(Kiso+1).intercept = b0_all;
ERR(Kiso+1).alpha = alpha_all;
ERR(Kiso+1).model = 'std(mu) = 1/(1+mu) - 1/(1+alpha*mu) + intercept  [global]';
ERR(Kiso+1).sigma_fun = sigma_fun_all;


%%% --- 2 --- Plot all six fits (5 isoforms + GLOBAL)
figure; 
tl = tiledlayout(3,2,'TileSpacing','compact'); 
for k = 1:Kiso+1
    nexttile; 
    if k <= Kiso
        x = mu_ex2(:,k); 
        y = std_ex2(:,k);
    else
        x = mu_ex2(:);   % GLOBAL uses all points
        y = std_ex2(:);
    end
    xfit = linspace(min(x), max(x), 200);
    yfit = ERR(k).sigma_fun(xfit);
    plot(x*100, y*100, 'o'); hold on;
    plot(xfit*100, yfit*100, 'k','LineWidth',1.3); grid on
    title(sprintf('%s | \\alpha=%.3g, b0=%.3g', ERR(k).name, ERR(k).alpha, ERR(k).intercept));
    xlabel('Mean (%)'); ylabel('Std (%)');
end

%%   K  with b0
load("datahekrep1nowt.mat");
load("datahekrep2nowt.mat");
load("datahekrep3nowt.mat");

datahekrep1 = datahekrep1nowt; 
datahekrep2 = datahekrep2nowt; 
datahekrep3 = datahekrep3nowt;



% Exon 1,2,3 
isoCols  = 2:6;                        % Inclusion, Skipping, FullIR, FirstIR, SecondIR
exonsUse = [2 3 4];                    % {2}=exon1, {3}=exon2, {4}=exon3

mu_all  = [];                          % (N1+N2+N3) x 5
std_all = [];

for ee = exonsUse
    [Cee, i1, i2, i3] = mintersect( ...
        datahekrep1{ee}(:,1), datahekrep2{ee}(:,1), datahekrep3{ee}(:,1));

    R1 = datahekrep1{ee}(i1, isoCols);   % Nie x 5
    R2 = datahekrep2{ee}(i2, isoCols);
    R3 = datahekrep3{ee}(i3, isoCols);

    X  = cat(3, R1, R2, R3);
    mu_e  = mean(X, 3);                  % Nie x 5
    std_e = std (X, 0, 3);               % Nie x 5

    mu_all  = [mu_all;  mu_e];
    std_all = [std_all; std_e];
end

mu_ex2  = mu_all;                    
std_ex2 = std_all; 

isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};
Kiso = numel(isoCols);

%%% --- 1 --- Fitting by isoform
ERR = struct('name',[],'alpha',[],'b0',[],'sigma_fun',[],'model',[]);
opts = optimoptions('lsqnonlin','Display','off');

for k = 1:Kiso
    x = mu_ex2(:,k);
    y = std_ex2(:,k);

    x = max(min(x, 1-1e-9), 1e-9);
    K = 1./x - 1;

    % Initialize intercept from linear fit
    p_lin = polyfit(x, y, 1);
    b0_init = max(p_lin(2), 0);
    alpha0 = 1;

    % Fit alpha & b0 jointly: σ(μ)=1/(1+αK) - 1/(1+K) + b0
    errfun = @(p) -y + (1./(1 + p(1).*K) - 1./(1 + K) + p(2));
    lb = [0, 0]; ub = [100, 1];
    p_fit = lsqnonlin(errfun, [alpha0, b0_init], lb, ub, opts);

    alpha_k = p_fit(1);
    b0_k = p_fit(2);

    % Define σ(μ) function
    sigma_fun_k = @(muq) max(1./(1 + alpha_k.*(1./muq - 1)) - 1./(1 + (1./muq - 1)) + b0_k, 0);

    ERR(k).name = isoNames{k};
    ERR(k).alpha = alpha_k;
    ERR(k).b0 = b0_k;
    ERR(k).model = 'σ(μ)=1/(1+αK)-1/(1+K)+b₀';
    ERR(k).sigma_fun = sigma_fun_k;
end

%%% --- 1A --- GLOBAL model
x_all = mu_ex2(:);
y_all = std_ex2(:);
x_all = max(min(x_all, 1-1e-9), 1e-9);
K_all = 1./x_all - 1;

p_lin_all = polyfit(x_all, y_all, 1);
b0_init_all = max(p_lin_all(2), 0);
alpha0 = 1;

errfun_all = @(p) -y_all + (1./(1 + p(1).*K_all) - 1./(1 + K_all) + p(2));
p_fit_all = lsqnonlin(errfun_all, [alpha0, b0_init_all], [0,0], [100,1], opts);

alpha_all = p_fit_all(1);
b0_all = p_fit_all(2);

sigma_fun_all = @(muq) max(1./(1 + alpha_all.*(1./muq - 1)) - 1./(1 + (1./muq - 1)) + b0_all, 0);

ERR(Kiso+1).name = 'GLOBAL';
ERR(Kiso+1).alpha = alpha_all;
ERR(Kiso+1).b0 = b0_all;
ERR(Kiso+1).model = 'σ(μ)=1/(1+αK)-1/(1+K)+b₀ [global]';
ERR(Kiso+1).sigma_fun = sigma_fun_all;

%%% --- 2 --- Plot all fits
% === Settings ===
Kiso = 5;
YMAX = 25;
isoColors = [
    0.0000 0.4470 0.7410;  % Inclusion (blue)
    0.8500 0.3250 0.0980;  % Skipping (orange)
    0.9290 0.6940 0.1250;  % FullIR (yellow)
    0.4940 0.1840 0.5560;  % FirstIR (purple)
    0.4660 0.6740 0.1880;  % SecondIR (green)
];

% === Layout size (in centimeters) ===
panel_w = 8.25;
panel_h = 8.5;
cols = 2;
rows = 3;
gap = 0.3;

fig_w = cols * panel_w + (cols - 1) * gap;
fig_h = rows * panel_h + (rows - 1) * gap;

fig = figure('Color','w');
set(fig, 'Units','centimeters', 'Position',[1 1 fig_w fig_h], ...
         'PaperUnits','centimeters', 'PaperPosition',[0 0 fig_w fig_h], ...
         'PaperSize',[fig_w fig_h]);

tl = tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');

for k = 1:Kiso+1
    ax = nexttile; hold on;

    if k <= Kiso
        x = mu_ex2(:,k);
        y = std_ex2(:,k);
        col = isoColors(k,:);
        name = ERR(k).name;
    else
        x = mu_ex2(:);
        y = std_ex2(:);
        col = [0 0 0];
        name = 'GLOBAL';
    end

    % scatter
    scatter(x*100, y*100, 20, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col, ...
        'MarkerFaceAlpha', 0.6);

    % fitted curve
    xfit = linspace(0,1,300);
    yfit = ERR(k).sigma_fun(xfit);
    plot(xfit*100, yfit*100, '-', ...
        'Color', col, 'LineWidth', 1.5);

    % axis formatting
    xlim([0 100]);
    ylim([0 YMAX]);
    box on;
    axis square;

    xlabel('Mean isoform frequency (%)', 'FontSize', 9);
    ylabel('Std (%)', 'FontSize', 9);

    title(sprintf('%s | \\alpha = %.3g, b_0 = %.3g', ...
        name, ERR(k).alpha, ERR(k).b0), ...
        'FontSize', 9, 'FontWeight','normal');
end

% === Save ===
fig_name = 'Figure_ErrorModel_sigma_vs_mu';
savefig(fig, [fig_name '.fig']);
print(fig, '-dsvg', '-r800', '-painters', [fig_name '.svg']);

print(gcf, '-dsvg', 'ErrorModel_sigma_vs_mu.svg');



%% save data
output_dir = getenv('ERROR_MODEL_OUTPUT_DIR');
if isempty(output_dir)
    output_dir = pwd;
end
outname_local = fullfile(output_dir, 'error_model_hek_withGLOBAL_allExons_withOffset.mat');

INFO = struct();
INFO.formula      = 'std(mu) = 1/(1+alpha*K) - 1/(1+K) + b0';
INFO.fit_space    = 'K-transform';
INFO.exonsUse     = exonsUse;
INFO.isoCols      = isoCols;
INFO.isoNames     = isoNames(:).';
INFO.hasIntercept = true;
INFO.date         = datestr(now, 31);
INFO.alpha_perIso = arrayfun(@(e)e.alpha, ERR(1:Kiso));
INFO.b0_perIso    = arrayfun(@(e)e.b0, ERR(1:Kiso));
INFO.alpha_global = ERR(Kiso+1).alpha;
INFO.b0_global    = ERR(Kiso+1).b0;
INFO.N_points     = size(mu_ex2,1);
INFO.note         = 'HEK293 3 replicates; exon1+2+3 intersect per exon; fitted with offset.';

save(outname_local, 'ERR', 'mu_ex2', 'std_ex2', 'isoNames', 'exonsUse', 'isoCols', 'INFO');

%%   K  with b0 -- all in one table (T,T1,T2,T3)
% Load replicates
load("datahekrep1nowt.mat");
load("datahekrep2nowt.mat");
load("datahekrep3nowt.mat");

datahekrep1 = datahekrep1nowt;
datahekrep2 = datahekrep2nowt;
datahekrep3 = datahekrep3nowt;

% Load updated error model (WITH offset)
load("error_model_hek_withGLOBAL_allExons_withOffset.mat");  % contains ERR(1:5) and ERR(6)

% Config
isoCols = 2:6;                          % Inclusion, Skipping, FullIR, FirstIR, SecondIR
isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};
exonsUse = [2 3 4];                     % Exon 1, 2, 3

% Init output tables
allRows = [];
T1 = table(); T2 = table(); T3 = table();

% --- Generate T (intersection of replicates) ---
for ee = exonsUse
    [Cee, i1, i2, i3] = mintersect( ...
        datahekrep1{ee}(:,1), datahekrep2{ee}(:,1), datahekrep3{ee}(:,1));

    R1 = datahekrep1{ee}(i1, isoCols);
    R2 = datahekrep2{ee}(i2, isoCols);
    R3 = datahekrep3{ee}(i3, isoCols);

    MU = (R1 + R2 + R3) / 3;

    STDspecific = nan(size(MU));
    for k = 1:numel(isoNames)
        STDspecific(:,k) = ERR(k).sigma_fun(MU(:,k));   % already includes offset b0
    end

    STDglobal = ERR(6).sigma_fun(MU);                    % global error model with offset

    tableBlock = [ ...
        Cee, ...
        R1, R2, R3, ...
        MU, ...
        STDspecific, ...
        STDglobal ...
    ];

    allRows = [allRows; tableBlock];
end

colNames = [{'Position'}, ...
            strcat('Rep1_', isoNames), ...
            strcat('Rep2_', isoNames), ...
            strcat('Rep3_', isoNames), ...
            strcat('Mean_', isoNames), ...
            strcat('StdSpec_', isoNames), ...
            strcat('StdGlob_', isoNames)];

T = array2table(allRows, 'VariableNames', colNames);


% --- Helper functions for replicate-specific tables ---
compute_STD_specific = @(mu_matrix) cell2mat(arrayfun( ...
    @(k) ERR(k).sigma_fun(mu_matrix(:,k)), ...
    1:numel(isoNames), 'UniformOutput', false));

compute_STD_global = @(mu_matrix) ERR(6).sigma_fun(mu_matrix);

% --- Generate T1, T2, T3 ---
for ee = exonsUse
    D1 = datahekrep1{ee};
    D2 = datahekrep2{ee};
    D3 = datahekrep3{ee};

    % Rep1
    mu1 = D1(:, isoCols);
    std1_spec = compute_STD_specific(mu1);
    std1_glob = compute_STD_global(mu1);
    tbl1 = array2table([D1(:,1), mu1, std1_spec, std1_glob], ...
        'VariableNames', [{'Position'}, ...
                          isoNames, ...
                          strcat('STDspec_', isoNames), ...
                          strcat('STDglob_', isoNames)]);
    T1 = [T1; tbl1];

    % Rep2
    mu2 = D2(:, isoCols);
    std2_spec = compute_STD_specific(mu2);
    std2_glob = compute_STD_global(mu2);
    tbl2 = array2table([D2(:,1), mu2, std2_spec, std2_glob], ...
        'VariableNames', [{'Position'}, ...
                          isoNames, ...
                          strcat('STDspec_', isoNames), ...
                          strcat('STDglob_', isoNames)]);
    T2 = [T2; tbl2];

    % Rep3
    mu3 = D3(:, isoCols);
    std3_spec = compute_STD_specific(mu3);
    std3_glob = compute_STD_global(mu3);
    tbl3 = array2table([D3(:,1), mu3, std3_spec, std3_glob], ...
        'VariableNames', [{'Position'}, ...
                          isoNames, ...
                          strcat('STDspec_', isoNames), ...
                          strcat('STDglob_', isoNames)]);
    T3 = [T3; tbl3];
end


writetable(T,  'D2D_input_T.csv');
writetable(T1, 'D2D_input_T1.csv');
writetable(T2, 'D2D_input_T2.csv');
writetable(T3, 'D2D_input_T3.csv');

%% .def & .csv(with offset)
% Load error model (with offset)
load('error_model_hek_withGLOBAL_allExons_withOffset.mat');  % has ERR(1:6)

% Load experimental PSI data
T = readtable('RON_HEK293.csv');  % existing file with mean values
mu = table2array(T(:, 2:6));      % Inclusion..SecondIR
isoNames = {'inclusion','skipping','fullir','firstir','secondir'};

% Compute STD from both error models
STDspec = nan(1, numel(isoNames));
STDglob = nan(1, numel(isoNames));

for k = 1:numel(isoNames)
    STDspec(k) = ERR(k).sigma_fun(mu(:,k));   % isoform-specific model
end
STDglob = ERR(6).sigma_fun(mu);               % global model

% Create readable table
result = table(isoNames', mu', STDspec', STDglob', ...
    'VariableNames', {'Isoform','Mean','STD_specific','STD_global'});
disp(result);

% === Option A: Write .def file (using specific STD) ===
def_file = 'RON_HEK293.def';
fid = fopen(def_file, 'w');

fprintf(fid, 'DESCRIPTION\n');
fprintf(fid, '"RON_HEK293_WT"\n\n');
fprintf(fid, 'PREDICTOR\n');
fprintf(fid, 't\tT\ts\ttime\t0\t100000000\n\n');
fprintf(fid, 'INPUTS\n\n');
fprintf(fid, 'OBSERVABLES\n\n');
fprintf(fid, 'ERRORS\n\n');
fprintf(fid, 'CONDITIONS\n');

for k = 1:numel(isoNames)
    fprintf(fid, 'sd_%s\t"%.6f"\n', isoNames{k}, STDspec(k));
end
fclose(fid);

fprintf('✅ DEF file written (specific σ): %s\n', def_file);

% === Option B: also write global-version (for comparison) ===
def_file_glob = 'RON_HEK293.def';
fid = fopen(def_file_glob, 'w');

fprintf(fid, 'DESCRIPTION\n');
fprintf(fid, '"RON_HEK293_WT"\n\n');
fprintf(fid, 'PREDICTOR\n');
fprintf(fid, 't\tT\ts\ttime\t0\t100000000\n\n');
fprintf(fid, 'INPUTS\n\n');
fprintf(fid, 'OBSERVABLES\n\n');
fprintf(fid, 'ERRORS\n\n');
fprintf(fid, 'CONDITIONS\n');

for k = 1:numel(isoNames)
    fprintf(fid, 'sd_%s\t"%.6f"\n', isoNames{k}, STDglob(k));
end
fclose(fid);

fprintf('✅ DEF file written (global σ): %s\n', def_file_glob);


%% replicate deviations from mean
load("datahekrep1.mat");
load("datahekrep2.mat");
load("datahekrep3.mat");

% Define helper functions
fix_one_column = @(x) round(floor(x) + (x - floor(x))*10, 4);
should_fix = @(x) (x - floor(x) > 0) & (x - floor(x) < 0.1 + 1e-12);

exonsUse = [2 3 4];

for ee = exonsUse
    pos = datahekrep2{ee}(:,1);
    mask = should_fix(pos);

    oldVals = pos(mask);
    newVals = fix_one_column(oldVals);
    pos(mask) = newVals;
    datahekrep2{ee}(:,1) = pos;
end

for ee = exonsUse
    pos = datahekrep3{ee}(:,1);
    mask = should_fix(pos);

    oldVals = pos(mask);
    newVals = fix_one_column(oldVals);
    pos(mask) = newVals;
    datahekrep3{ee}(:,1) = pos;
end

% Exon 1,2,3 
isoCols  = 2:6;                        % Inclusion, Skipping, FullIR, FirstIR, SecondIR
isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};
exonsUse = [2 3 4];                    % {2}=exon1, {3}=exon2, {4}=exon3

allRows = [];

for ee = exonsUse
    [Cee, i1, i2, i3] = mintersect( ...
        datahekrep1{ee}(:,1), datahekrep2{ee}(:,1), datahekrep3{ee}(:,1));

    R1 = datahekrep1{ee}(i1, isoCols);   % Nie x 5
    R2 = datahekrep2{ee}(i2, isoCols);
    R3 = datahekrep3{ee}(i3, isoCols);

   % --- Compute mean values ---
    MU = (R1 + R2 + R3) / 3;

    % --- (optional) if you already have error model, compute std ---
    STD = std(cat(3, R1, R2, R3), 0, 3);

    % --- Combine all columns into one block ---
    tableBlock = [ ...
        Cee, ...             % position
        R1, R2, R3, ...      % replicates
        MU, STD ...          % mean + std
    ];

    % Append to total list
    allRows = [allRows; tableBlock];
end

% --- Define column names ---
colNames = [{'Position'}, ...
            strcat('Rep1_', isoNames), ...
            strcat('Rep2_', isoNames), ...
            strcat('Rep3_', isoNames), ...
            strcat('Mean_', isoNames), ...
            strcat('Std_', isoNames)];

% --- Convert to table ---
T = array2table(allRows, 'VariableNames', colNames);

% --- Save as CSV and MAT file ---
writetable(T, 'T.csv');



% ============= Load table  ============= 
T = readtable('T.csv');

% Define isoform names
isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};

% Prepare figure folders
outDir = 'Isoform_Deviation_Plots';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ===  isoform  ===
for k = 1:numel(isoNames)
    iso = isoNames{k};

    %  mean & replicate 
    mean_col = T.(sprintf('Mean_%s', iso));
    rep1_col = T.(sprintf('Rep1_%s', iso));
    rep2_col = T.(sprintf('Rep2_%s', iso));
    rep3_col = T.(sprintf('Rep3_%s', iso));

    % Δ = Rep - Mean
    delta1 = rep1_col - mean_col;
    delta2 = rep2_col - mean_col;
    delta3 = rep3_col - mean_col;

    % plot
    figure('Name', iso, 'Color', 'w');
    hold on;
    scatter(mean_col, delta1, 40, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(mean_col, delta2, 40, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(mean_col, delta3, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.6);

    % lines
    yline(0, '--k', 'LineWidth', 1);
    xlabel(sprintf('Mean %s', iso));
    ylabel('Replicate − Mean');
    title(sprintf('Replicate deviations for %s', iso));
    legend({'Rep1','Rep2','Rep3'}, 'Location', 'best');
    grid on;

    % save
    saveas(gcf, fullfile(outDir, sprintf('%s_deviation.png', iso)));

    savefig(gcf, fullfile(outDir, sprintf('%s_deviation.fig', iso)));
    print(gcf, fullfile(outDir, sprintf('%s_deviation.svg', iso)), '-dsvg');
    close;
end



% === isoform overall ===
figure('Name', 'All Isoforms', 'Color', 'w');
hold on;

colors = lines(numel(isoNames));
for k = 1:numel(isoNames)
    iso = isoNames{k};
    mean_col = T.(sprintf('Mean_%s', iso));
    rep1_col = T.(sprintf('Rep1_%s', iso));
    rep2_col = T.(sprintf('Rep2_%s', iso));
    rep3_col = T.(sprintf('Rep3_%s', iso));

    delta_all = [ ...
        rep1_col - mean_col; ...
        rep2_col - mean_col; ...
        rep3_col - mean_col];
    mean_all = [ ...
        mean_col; ...
        mean_col; ...
        mean_col];

    scatter(mean_all, delta_all, 25, 'filled', ...
        'MarkerFaceColor', colors(k,:), 'MarkerFaceAlpha', 0.5, ...
        'DisplayName', iso);
end

yline(0, '--k', 'LineWidth', 1);
xlabel('Mean PSI');
ylabel('Replicate − Mean');
title('Replicate deviations across all isoforms');
legend('Location','best');


saveas(gcf, fullfile(outDir, 'AllIsoforms_deviation.png'));

savefig(gcf, fullfile(outDir, sprintf('AllIsoforms_deviation.fig', iso)));
print(gcf, fullfile(outDir, sprintf('AllIsoforms_deviation.svg', iso)), '-dsvg');




%% Z deviation plots (log-ratio metric)
outDirZ = 'Isoform_LogRatio_Plots';
if ~exist(outDirZ, 'dir')
    mkdir(outDirZ);
end

isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};
nIso = numel(isoNames);

% Preallocate
Zrep1 = zeros(height(T), nIso);
Zrep2 = zeros(height(T), nIso);
Zrep3 = zeros(height(T), nIso);
Zmean = zeros(height(T), nIso);

% Small epsilon to avoid log10(0)
eps_val = 1e-9;

% --- Compute log10 ratios relative to Inclusion ---
for k = 1:nIso
    iso = isoNames{k};

    rep1 = T.(sprintf('Rep1_%s', iso));
    rep2 = T.(sprintf('Rep2_%s', iso));
    rep3 = T.(sprintf('Rep3_%s', iso));

    inc1 = T.Rep1_Inclusion;
    inc2 = T.Rep2_Inclusion;
    inc3 = T.Rep3_Inclusion;

    % log10 ratio to Inclusion
    % Zrep1(:,k) = log10( max(rep1 ./ max(inc1, eps_val), eps_val) );
    % Zrep2(:,k) = log10( max(rep2 ./ max(inc2, eps_val), eps_val) );
    % Zrep3(:,k) = log10( max(rep3 ./ max(inc3, eps_val), eps_val) );

    Zrep1(:,k) = log10( rep1 ./ inc1);
    Zrep2(:,k) = log10( rep2 ./ inc2);
    Zrep3(:,k) = log10( rep3 ./ inc3);

    % Mean across replicates
    Zmean(:,k) = mean([Zrep1(:,k), Zrep2(:,k), Zrep3(:,k)], 2);
end


% ---  Per-isoform deviation plots ---
for k = 1:nIso
    iso = isoNames{k};
    Zm = Zmean(:,k);
    d1 = Zrep1(:,k) - Zm;
    d2 = Zrep2(:,k) - Zm;
    d3 = Zrep3(:,k) - Zm;

    figure('Name', ['LogRatio_', iso], 'Color', 'w');
    hold on;
    scatter(Zm, d1, 40, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(Zm, d2, 40, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(Zm, d3, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.6);

    yline(0, '--k', 'LineWidth', 1);
    xlabel(sprintf('Mean log10(%s / Inclusion)', iso));
    ylabel('Replicate − Mean');
    title(sprintf('Deviation of log10(%s / Inclusion)', iso));
    legend({'Rep1','Rep2','Rep3'}, 'Location', 'best');
    grid on;

    saveas(gcf, fullfile(outDirZ, sprintf('%s_logRatio_dev.png', iso)));

    savefig(gcf, fullfile(outDirZ, sprintf('%s_logRatio_dev.fig', iso)));
    print(gcf, fullfile(outDirZ, sprintf('%s_logRatio_dev.svg', iso)), '-dsvg');
    close;
end


% ---  Combined deviation plot (all isoforms) ---
figure('Name', 'All Isoforms Log10 Ratios', 'Color', 'w');
hold on;
colors = lines(nIso);

for k = 1:nIso
    iso = isoNames{k};
    Zm = Zmean(:,k);
    d_all = [Zrep1(:,k)-Zm; Zrep2(:,k)-Zm; Zrep3(:,k)-Zm];
    Z_all = [Zm; Zm; Zm];

    scatter(Z_all, d_all, 25, 'filled', ...
        'MarkerFaceColor', colors(k,:), 'MarkerFaceAlpha', 0.5, ...
        'DisplayName', iso);
end

yline(0, '--k', 'LineWidth', 1);
xlabel('Mean log10(isoform / Inclusion)');
ylabel('Replicate − Mean');
title('Deviation of log10 ratios relative to Inclusion');
legend('Location','best');
grid on;

saveas(gcf, fullfile(outDirZ, 'AllIsoforms_LogRatio_deviation.png'));

savefig(gcf, fullfile(outDirZ, sprintf('AllIsoforms_LogRatio_deviation.fig', iso)));
print(gcf, fullfile(outDirZ, sprintf('AllIsoforms_LogRatio_deviation.svg', iso)), '-dsvg');


%% Check if 5 isoforms sum to 1 (Rep1/Rep2/Rep3/Mean)
T = readtable('D2D_input_T.csv');

iso = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};
eps_tol = 1e-6;   % tolerance for "close to 1"

% Helper to compute row-wise sums and deviations
rowSum = @(prefix) sum( [ T.(sprintf('%s_%s',prefix,iso{1})) , ...
                          T.(sprintf('%s_%s',prefix,iso{2})) , ...
                          T.(sprintf('%s_%s',prefix,iso{3})) , ...
                          T.(sprintf('%s_%s',prefix,iso{4})) , ...
                          T.(sprintf('%s_%s',prefix,iso{5})) ] , 2);

S1 = rowSum('Rep1');
S2 = rowSum('Rep2');
S3 = rowSum('Rep3');
SM = rowSum('Mean');

load("datahekrep1nowt.mat");
datahekrep1 = datahekrep1nowt;

isoCols = 2:6;          % Inclusion, Skipping, FullIR, FirstIR, SecondIR
exonsUse = [2 3 4];     

for ee = exonsUse
    exonData = datahekrep1{ee};
    sumIso = sum(exonData(:, isoCols), 2);
    
    % 
    exonData(:, 7) = sumIso;

    % 
    datahekrep1{ee} = exonData; 
end



%% Replicate deviation analysis with Gaussian fits
% Load data
T = readtable('T.csv');
isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};

% Prepare output folder
outDir = 'Isoform_Deviation_Plots';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% --- Helper function for computing deviations and Gaussian fit ---
getDeviations = @(T, iso) [ ...
    T.(sprintf('Rep1_%s', iso)) - T.(sprintf('Mean_%s', iso)); ...
    T.(sprintf('Rep2_%s', iso)) - T.(sprintf('Mean_%s', iso)); ...
    T.(sprintf('Rep3_%s', iso)) - T.(sprintf('Mean_%s', iso)) ];

% --- Individual plots per isoform ---
for k = 1:numel(isoNames)
    iso = isoNames{k};
    delta_all = getDeviations(T, iso);
    mu = mean(delta_all);
    sigma = std(delta_all);

    % Histogram (normalized to PDF)
    figure('Name', iso, 'Color', 'w'); hold on;
    histogram(delta_all, 40, 'Normalization', 'pdf', ...
        'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    % Gaussian fit
    xfit = linspace(min(delta_all), max(delta_all), 200);
    yfit = normpdf(xfit, mu, sigma);
    plot(xfit, yfit, 'r-', 'LineWidth', 2);

    % Formatting
    xline(0, '--k', 'LineWidth', 1);
    xlabel('Replicate − Mean');
    ylabel('Density');
    title(sprintf('%s deviation (μ = %.3g, σ = %.3g)', iso, mu, sigma));
    legend({'Histogram','Gaussian fit'}, 'Location', 'best');
    grid on;

    saveas(gcf, fullfile(outDir, sprintf('%s_deviation_fit.png', iso)));

    savefig(gcf, fullfile(outDir, sprintf('%s_deviation_fit.fig', iso)));
    print(gcf, fullfile(outDir, sprintf('%s_deviation_fit.svg', iso)), '-dsvg');
    close;
end

% --- Combined Gaussian fit curves ---
fprintf('=== Generating combined Gaussian curve comparison ===\n');
colors = lines(numel(isoNames));
figure('Name', 'All Isoforms Gaussian Fits', 'Color', 'w'); hold on;

for k = 1:numel(isoNames)
    iso = isoNames{k};
    delta_all = getDeviations(T, iso);
    mu = mean(delta_all);
    sigma = std(delta_all);
    xfit = linspace(-0.5, 0.5, 400);
    yfit = normpdf(xfit, mu, sigma);

    plot(xfit, yfit, 'LineWidth', 2, 'Color', colors(k,:), ...
        'DisplayName', sprintf('%s (σ=%.3g)', iso, sigma));
end

xline(0, '--k', 'LineWidth', 1);
xlabel('Replicate − Mean');
ylabel('Density');
title('Gaussian-fitted deviation distributions across isoforms');
legend('Location', 'best'); grid on;

saveas(gcf, fullfile(outDir, 'AllIsoforms_GaussianFits.png'));

savefig(gcf, fullfile(outDir, sprintf('AllIsoforms_GaussianFits.fig', iso)));
print(gcf, fullfile(outDir, sprintf('AllIsoforms_GaussianFits.svg', iso)), '-dsvg');

% --- Combined histogram + Gaussian fit overlay ---
fprintf('=== Generating combined histogram + Gaussian fit comparison ===\n');
figure('Name', 'All Isoforms Gaussian Fits + Histograms', 'Color', 'w'); hold on;

for k = 1:numel(isoNames)
    iso = isoNames{k};
    delta_all = getDeviations(T, iso);
    mu = mean(delta_all);
    sigma = std(delta_all);

    % Histogram (semi-transparent)
    histogram(delta_all, 40, 'Normalization', 'pdf', ...
        'FaceColor', colors(k,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
        'DisplayName', sprintf('%s data', iso));

    % Gaussian fit
    xfit = linspace(min(delta_all), max(delta_all), 400);
    yfit = normpdf(xfit, mu, sigma);
    plot(xfit, yfit, '-', 'Color', colors(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s fit (σ=%.3g)', iso, sigma));

    % Optional: mark mean
    plot(mu, normpdf(mu, mu, sigma), 'o', 'Color', colors(k,:), ...
        'MarkerFaceColor', colors(k,:), 'MarkerSize', 5, 'HandleVisibility', 'off');
end

xline(0, '--k', 'LineWidth', 1);
xlabel('Replicate − Mean');
ylabel('Density');
title('Deviation distributions across isoforms (histogram + Gaussian fit)');
legend('Location', 'best'); grid on;

saveas(gcf, fullfile(outDir, 'AllIsoforms_GaussianFits_withHist.png'));

savefig(gcf, fullfile(outDir, sprintf('AllIsoforms_GaussianFits_withHist.fig', iso)));
print(gcf, fullfile(outDir, sprintf('AllIsoforms_GaussianFits_withHist.svg', iso)), '-dsvg');

%% D2D_input_T1 with WT
% --- Load T1 (without WT) ---
T1 = readtable('D2D_input_T1.csv');

% --- Load WT source + error model ---
load("datahekrep1nowt.mat");  % WT is stored in datahekrep1{1}
load("error_model_hek_withGLOBAL_allExons_withOffset.mat");  % contains ERR

% --- Settings ---
isoCols  = 2:6;  % Inclusion, Skipping, FullIR, FirstIR, SecondIR
isoNames = {'Inclusion','Skipping','FullIR','FirstIR','SecondIR'};

% --- Helper function ---
compute_STD_specific = @(mu_matrix) cell2mat(arrayfun( ...
    @(k) ERR(k).sigma_fun(mu_matrix(:,k)), ...
    1:numel(isoNames), 'UniformOutput', false));

% ---------------------- WT ----------------------
D1_WT = datahekrep1nowt{1};                 % WT data (exon0)
pos_wt = D1_WT(:,1);                        % Position (1 x 1)
mu_wt = D1_WT(:, isoCols);                 % 1 x 5
std_spec_wt = compute_STD_specific(mu_wt); % 1 x 5
row_wt = [pos_wt, std_spec_wt];            % 1 x 6

% ---------------------- KD ----------------------
% Extract [Position, STDspec_*] columns from T1
std_cols = strcat('STDspec_', isoNames);
T1_extract = T1(:, ['Position', std_cols]);
rows_kd = table2array(T1_extract);         % 320 x 6

% ---------------------- Concatenate & Save ----------------------
error_fitpred_new = [row_wt; rows_kd];     % 321 x 6 double
save('error_fitpred_new.mat', 'error_fitpred_new');

%% error model -- alphas

mu = linspace(1e-4, 1-1e-4, 1000);
K  = (1./mu) - 1;

b0 = 0;
alphas = [0.1 0.8 1 1.25 10];

figure; hold on;

for a = alphas
    sigma = 1./(1 + a*K) - 1./(1 + K) + b0;

    plot(mu*100, sigma*100, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\alpha=%.2g', a));
end

xlabel('\mu (%)');
ylabel('\sigma(\mu) (%)');
title('Effect of \alpha on error model shape', 'FontWeight','normal');
legend('Location','best');

xticks(0:20:100);
yticks(-60:20:60);
box on;
xlim([0 100]);
ylim([-60 60]);






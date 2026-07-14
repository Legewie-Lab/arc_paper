%% analyse_fittingV2.m
% Post-processing and plotting for fitted one-step and two-step RON D2D models.
% Requires the fitted ar objects and replicate data files to be present in the
% MATLAB working directory.

%% Load data
load('hekrep1nowt.mat');
load('datahekrep1.mat');
% arLoad('Results/20251005T135917_two_step_2000_parallelFit');
load('multistart1_500_05.mat');
ar_one = ar;  
load('multistart2_500_05.mat');
ar_two = best_ar;

lenk1 = 62;
lenk2 = 274;
lenk3 = 69;
isoNames = {'inclusion','skipping','fullir','firstir','secondir'};

% Scale datar2 values to percentage
datar2(:, 2:6) = datar2(:, 2:6) .* 100;

% === Plot 1: Model vs Data (two panels: one-step and two-step) ===
figure('Name','Model_vs_Data_Comparison','Color','w');

models = {'One-step','Two-step'};
ars = {ar_one, ar_two};

for m = 1:2
    ar = ars{m};
    dr = ar.model.data;
    dt = {dr.yExp};
    prd = {dr.yExpSimu};
    titre = ["Inclusion", "Skipping", "FullIR", "1IR", "2IR"];

    % collect
    nData = length(dt);
    data = zeros(nData, size(dt{1}, 2));
    pred = zeros(nData, size(prd{1}, 2));
    for i = 1:nData
        data(i, :) = dt{i};
        pred(i, :) = prd{i};
    end

    % === Subplot ===
    subplot(1,2,m); hold on;
    
    % --- scatter: Model vs Data ---
    h1 = plot(data .* 100, pred .* 100, 'o', 'MarkerSize',4);
    
    isoNames = {'Inclusion','Skipping','FullIR','1IR','2IR'};
    for i = 1:min(length(h1), numel(isoNames))
        h1(i).DisplayName = isoNames{i};
    end
    
    % --- background uncertainty band (no legend) ---
    p  = [-0.0894, 0.0914, 0.0140];
    yw = linspace(0,1,200);
    vv = polyval(p, yw);
    
    fill([yw fliplr(yw)]*100, ...
         [yw-vv fliplr(yw+vv)]*100, ...
         [0.7 0.7 0.7], ...
         'FaceAlpha',0.4, ...
         'EdgeColor','none', ...
         'HandleVisibility','off');
    
    % --- diagonal ---
    hdiag = plot([0 100],[0 100],'k-','LineWidth',1.25);
    hdiag.DisplayName = 'Perfect fit';
    
    % --- colors ---
    colorOrder = get(gca,'ColorOrder');
    for i = 1:min(size(data,2), size(pred,2))
        set(h1(i),'Color',colorOrder(i,:));
    end
    
    % --- axes & legend ---
    xlabel('Data (%)');
    ylabel('Model (%)');
    xlim([0 100]); ylim([0 100]);
    axis square;box on; 
    title(models{m},'FontWeight','bold');
    legend('Location','southeast');
end

% exportgraphics(gcf, ...
%     'Plot1_Model_vs_Data_Both.png', ...
%     'Resolution', 300);

print(gcf, '-dsvg', 'Plot1_Model_vs_Data_Both.svg');
savefig(gcf, 'Plot1_Model_vs_Data_Both.fig');


% === Plot 1.1 : Model vs Data (exon 2 only & exon 1 3) ===
figure('Name','Model_vs_Data_Exon2','Color','w');
for m = 1:2
    ar = ars{m};
    dr = ar.model.data;
    dt = {dr.yExp};
    prd = {dr.yExpSimu};

    % exon2 
    exon2_idx = (lenk1+2):(lenk1+lenk2+1);

    % collect
    data = zeros(length(exon2_idx), size(dt{1}, 2));
    pred = zeros(length(exon2_idx), size(prd{1}, 2));
    for ii = 1:length(exon2_idx)
        data(ii, :) = dt{exon2_idx(ii)};
        pred(ii, :) = prd{exon2_idx(ii)};
    end

    subplot(1,2,m); hold on;
    h1 = plot(data .* 100, pred .* 100, 'o', 'MarkerSize', 4);

    % 
    p  = [-0.0894, 0.0914, 0.0140];
    yw = linspace(0,1,200);
    vv = polyval(p, yw);
    fill([yw fliplr(yw)]*100, [yw-vv fliplr(yw+vv)]*100, ...
        [0.7 0.7 0.7], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'none', 'HandleVisibility','off');

    plot([0 100],[0 100],'k-','LineWidth',1.25);

    colorOrder = get(gca,'ColorOrder');
    for i = 1:min(size(data,2), size(pred,2))
        set(h1(i),'Color',colorOrder(i,:));
    end

    xlabel('Data (%)'); ylabel('Model (%)');
    xlim([0 100]); ylim([0 100]); axis square; box on;
    title(models{m} + " (Exon2 only)",'FontWeight','bold');
end

print(gcf, '-dsvg', 'Plot1.1_Model_vs_Data_Exon2.svg');
savefig(gcf, 'Plot1.1_Model_vs_Data_Exon2.fig');



figure('Name','Model_vs_Data_Exon1_3','Color','w');
for m = 1:2
    ar = ars{m};
    dr = ar.model.data;
    dt = {dr.yExp};
    prd = {dr.yExpSimu};

    exon1_idx = 2:(1+lenk1);
    exon3_idx = (lenk1+lenk2+2):(lenk1+lenk2+lenk3+1);
    selected_idx = [exon1_idx, exon3_idx];

    data = zeros(length(selected_idx), size(dt{1}, 2));
    pred = zeros(length(selected_idx), size(prd{1}, 2));
    for ii = 1:length(selected_idx)
        data(ii, :) = dt{selected_idx(ii)};
        pred(ii, :) = prd{selected_idx(ii)};
    end

    subplot(1,2,m); hold on;
    h1 = plot(data .* 100, pred .* 100, 'o', 'MarkerSize', 4);

    p  = [-0.0894, 0.0914, 0.0140];
    yw = linspace(0,1,200);
    vv = polyval(p, yw);
    fill([yw fliplr(yw)]*100, [yw-vv fliplr(yw+vv)]*100, ...
        [0.7 0.7 0.7], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'none', 'HandleVisibility','off');

    plot([0 100],[0 100],'k-','LineWidth',1.25);

    colorOrder = get(gca,'ColorOrder');
    for i = 1:min(size(data,2), size(pred,2))
        set(h1(i),'Color',colorOrder(i,:));
    end

    xlabel('Data (%)'); ylabel('Model (%)');
    xlim([0 100]); ylim([0 100]); axis square; box on;
    title(models{m} + " (Exon1 & 3)",'FontWeight','bold');
end

print(gcf, '-dsvg', 'Plot1.2_Model_vs_Data_Exon1_3.svg');
savefig(gcf, 'Plot1.2_Model_vs_Data_Exon1_3.fig');


% === Plot 2.1: PSI vs Efficiency comparison (Data + One-step + Two-step) ===
figure('Color','w'); hold on;

% --- Experimental data ---
psi_mut = datahekrep1{3}(:,2) ./ (datahekrep1{3}(:,2) + datahekrep1{3}(:,3)) * 100;
eff_mut = (datahekrep1{3}(:,2) + datahekrep1{3}(:,3)) * 100;
psi_wt  = datahekrep1{1}(:,2) ./ (datahekrep1{1}(:,2) + datahekrep1{1}(:,3)) * 100;
eff_wt  = (datahekrep1{1}(:,2) + datahekrep1{1}(:,3)) * 100;

plot(psi_mut, eff_mut, 'go', 'DisplayName','Mutation data');
plot(psi_wt,  eff_wt,  'k.', 'MarkerSize',25, 'DisplayName','WT data');

% --- index of k2 region (same as your old valk2) ---
idx_k2 = lenk1 + 2 : lenk1 + lenk2 + 1;

% --- Helper: extract predictions safely ---
extract_pred = @(arS) vertcat(arS.model.data(:).yExpSimu);   % works for both scalar or array structs

% --- One-step model ---
P_one   = extract_pred(ar_one);
valk2_1 = P_one(idx_k2, :);
da_one  = [ valk2_1(:,1)./(valk2_1(:,1)+valk2_1(:,2)), 1 - sum(valk2_1(:,3:5),2) ] * 100;
g_one   = sortrows(da_one,1);
% Scatter plot
scatter(da_one(:,1), da_one(:,2), 25, 'b', 'filled', 'MarkerFaceAlpha',0.3, 'DisplayName','One-step points');
% Fitting curve
plot(g_one(:,1), g_one(:,2), 'b-', 'LineWidth',2, 'DisplayName','One-step model');

% --- Two-step model ---
P_two   = extract_pred(ar_two);
valk2_2 = P_two(idx_k2, :);
da_two  = [ valk2_2(:,1)./(valk2_2(:,1)+valk2_2(:,2)), 1 - sum(valk2_2(:,3:5),2) ] * 100;
g_two   = sortrows(da_two,1);
% Scatter plot
scatter(da_two(:,1), da_two(:,2), 25, 'r', 'filled', 'MarkerFaceAlpha',0.3, 'DisplayName','Two-step points');
% Fitting curve
plot(g_two(:,1), g_two(:,2), 'r-', 'LineWidth',2, 'DisplayName','Two-step model');

% --- Aesthetics ---
xlabel('PSI (%)'); ylabel('Efficiency (%)');
xlim([0 100]); ylim([40 100]);
title('PSI–Efficiency: Data vs One-step vs Two-step');
legend('Location','bestoutside'); axis square; 

% exportgraphics(gcf, ...
%     'Plot2_PSI_vs_Efficiency_exon2.png', ...
%     'Resolution', 300);

print(gcf, '-dsvg', 'Plot2_PSI_vs_Efficiency_exon2.svg');
savefig(gcf, 'Plot2_PSI_vs_Efficiency_exon2.fig');




%% === Chi² analysis for both models (from D2D & error model, including k1/k2/k3 groups) ===
models = {'One-step','Two-step'};
ars = {ar_one, ar_two};

for m = 1:2
    ar = ars{m};
    fprintf('\n==================== %s ====================\n', models{m});

    % === collect model outputs ===
    dr = ar.model.data;
    chis = zeros(length(dr), 5);
    data = zeros(length(dr), size(dr(1).yExp,2));
    pred = zeros(length(dr), size(dr(1).yExpSimu,2));
    for i = 1:length(dr)
        chis(i,:) = dr(i).chi2;
        data(i,:) = dr(i).yExp;
        pred(i,:) = dr(i).yExpSimu;
    end

    % --- chi² from D2D (global) ---
    chisum = sum(chis(:));           
    chi_mean = mean(chis(:));    
    fprintf('Total chi2 (D2D internal) = %.2f\n', chisum);
    fprintf('Mean chi2  = %.2f\n', chi_mean);

    % === group indices ===
    idx_k1 = 2:lenk1+1;
    idx_k2 = lenk1+2 : lenk1+lenk2+1;
    idx_k3 = size(data,1) - lenk3 + 1 : size(data,1);
    chi_k1 = sum(chis(idx_k1, :), 'all');
    chi_k2 = sum(chis(idx_k2, :), 'all');
    chi_k3 = sum(chis(idx_k3, :), 'all');
    fprintf('Chi2 for k1 mutations = %.2f\n', chi_k1);
    fprintf('Chi2 for k2 mutations = %.2f\n', chi_k2);
    fprintf('Chi2 for k3 mutations = %.2f\n', chi_k3);

    % --- chi² from error model ---
    T1 = readtable('D2D_input_T1.csv');

    std_cols = {'STDspec_Inclusion','STDspec_Skipping','STDspec_FullIR','STDspec_FirstIR','STDspec_SecondIR'};
    % std_cols = {'STDglob_Inclusion','STDglob_Skipping','STDglob_FullIR','STDglob_FirstIR','STDglob_SecondIR'};
    
    data_mut = data(2:end, :);
    pred_mut = pred(2:end, :);
    std_matrix = T1{:, std_cols};
    N = size(data_mut, 1);
    chi2_all = zeros(N, 1);

    for i = 1:N
        obs = data_mut(i,:);
        exp = pred_mut(i,:);
        stdv = std_matrix(i,:);
        stdv(stdv==0) = 1e-6;
        chi2_all(i) = sum(((obs - exp)./stdv).^2);
    end

    % --- overall ---
    num_params = sum(ar.qFit == 1);
    num_points = N * size(data_mut,2);
    dof_total = num_points - num_params;
    chi2_total = sum(chi2_all);
    p_total = 1 - chi2cdf(chi2_total, dof_total);
    reduced_total = chi2_total / dof_total;

    fprintf('\n=== Overall Chi² (error model) ===\n');
    fprintf('  Total χ² = %.2f\n', chi2_total);
    fprintf('  Reduced χ² = %.2f\n', reduced_total);
    fprintf('  p-value = %.3e\n', p_total);

    % --- per group ---
    idx1 = 1:lenk1;
    idx2 = lenk1+1 : lenk1+lenk2;
    idx3 = size(data_mut,1) - lenk3 + 1 : size(data_mut,1);

    group_chi2 = [sum(chi2_all(idx1)); sum(chi2_all(idx2)); sum(chi2_all(idx3))];
    group_n = [numel(idx1); numel(idx2); numel(idx3)];
    group_data_points = group_n * size(data_mut, 2);
    group_params_used = [lenk1; lenk2; lenk3];
    group_dof = group_data_points - group_params_used;
    group_p = 1 - chi2cdf(group_chi2, group_dof);
    group_reduced = group_chi2 ./ group_dof;
    group_names = {'k1','k2','k3'};

    fprintf('\n=== Group-wise Chi² summary ===\n');
    for g = 1:3
        fprintf('\nGroup: %s\n', group_names{g});
        fprintf('  Data points       : %d\n', group_data_points(g));
        fprintf('  Fitted parameters : %d\n', group_params_used(g));
        fprintf('  Degrees of freedom: %d\n', group_dof(g));
        fprintf('  Total chi²        : %.2f\n', group_chi2(g));
        fprintf('  Reduced chi²      : %.2f\n', group_reduced(g));
        fprintf('  P-value           : %.4e\n', group_p(g));
    end

    % --- Save results table ---
    Chi2Summary = table( ...
        [group_names'; "All"], ...
        [group_chi2; chi2_total], ...
        [group_reduced; reduced_total], ...
        [group_p; p_total], ...
        'VariableNames', {'Group','Chi2','Reduced_Chi2','P_Value'});

    fname = sprintf('Chi2_summary_%s.csv', models{m});
    writetable(Chi2Summary, fname);
end










% Load updated model result
load('Controlnowt.mat');         
load('control1nowt_datar2.mat');        
load('datakontrolle.mat');

% Scale datar2 values to percentage
datar2(:, 2:6) = datar2(:, 2:6) .* 100;

% Define mutation counts
lenk1 = 115; lenk2 = 320; lenk3 = 139;

% Extract model data
dr = ar.model.data;
dt = {dr.yExp};
prd = {dr.yExpSimu};
ksolu = {dr.pNum};
at = array2table([dr(1).pNum]);
at.Properties.VariableNames = [dr(1).p];

% Parameter reordering
emr = {'k1', 'k2a', 'k2b', 'k3', 'kret', 'k4', 'k5a', 'k5b', 'k6', 'kspli', 'kincl', 'kskip', 'kdr1', 'kdr2', 's'};
aa = at(:, emr);  

% Initialize data collection
nData = length(dt);
data = zeros(nData, size(dt{1}, 2));
pred = zeros(nData, size(prd{1}, 2));
chis = zeros(nData, 5);

for i = 1:nData
    data(i, :) = dt{i};
    pred(i, :) = prd{i};
    chis(i,:) = dr(i).chi2;
end

chisum = sum(chis);

% Parameter extraction
k1_res = zeros(lenk1, 1);
numk1 = zeros(lenk1, 1);
k2_res = zeros(lenk2, 1);
numk2 = zeros(lenk2, 1);
k3_res = zeros(lenk3, 1);
numk3 = zeros(lenk3, 1);

for i = 1:lenk1
    k1_res(i) = ksolu{i + 1}(1);
    numk1(i) = datar2(i + 1, 1);
end

for i = 1:lenk2
    idx = i + lenk1 + 1;
    k2_res(i) = ksolu{idx}(1);
    numk2(i) = datar2(idx, 1);
end

for i = 1:lenk3
    idx = size(datar2, 1) - lenk3 + i;
    k3_res(i) = ksolu{idx}(1);
    numk3(i) = datar2(idx, 1);
end

%%%%%%% for k2_res %%%%%%%
% 1) k2_res: prepend WT = 1
k2_res_new = [1; k2_res];    % 321 x 1
k2_res = k2_res_new;

% 2) pred: keep WT + exon2 only
wt_idx   = 1;
k2_start = lenk1 + 2;              % 117
k2_end   = lenk1 + lenk2 + 1;      % 436
pred_new = pred([wt_idx, k2_start:k2_end], :);% 321 x 5
pred = pred_new;

% 3) aa: WT kinetic parameters (already correct)

% 4) Save as ONE file (drop-in replacement)
save('k2_res_new.mat', 'k2_res', 'pred', 'aa');

save('aa.mat', 'aa');


%% Control & KD Version 2
% === Load fitted model results and required variables ===
load('Controlnowt.mat');              % contains ar and ar.model.data
load('control1nowt_datar2.mat');      % contains datar2
load('datakontrolle.mat');            % optional for numk2 if needed

% scale datar2 values
datar2(:, 2:6) = datar2(:, 2:6) .* 100;

% counts of mutations
lenk1 = 115;
lenk2 = 320;
lenk3 = 139;

% === Extract WT kinetic parameters ===
dr = ar.model.data;
ksolu = {dr.pNum};
at = array2table([dr(1).pNum]);
at.Properties.VariableNames = [dr(1).p];
param_names = {'k1', 'k2a', 'k2b', 'k3', 'kret', 'k4', 'k5a', 'k5b', ...
               'k6', 'kspli', 'kincl', 'kskip', 'kdr1', 'kdr2', 's'};
aa = at(:, param_names);  % WT kinetic parameters

% === Extract f_mutX parameter values ===
f_mut_indices = find(contains(ar.pLabel, 'f_mut'));

labels = ar.pLabel(f_mut_indices);
mutation_nums = zeros(length(labels),1);
for i = 1:length(labels)
    str = labels{i};
    num_str = regexp(str, '\d+', 'match');
    if ~isempty(num_str)
        mutation_nums(i) = str2double(num_str{1});
    end
end

[~, sorted_idx] = sort(mutation_nums);
sorted_indices = f_mut_indices(sorted_idx);  % correctly sorted

k2_res = zeros(length(sorted_indices), 1);
for i = 1:length(sorted_indices)
    idx = sorted_indices(i);
    k2_res(i) = ksolu{idx}(1);   % get fitted scaling value for k2a/k2b
end
k2_res = [1; k2_res];  % prepend WT value = 1

% === Extract predictions in the same order ===
sorted_labels = ar.pLabel(sorted_indices);  % f_mut1 ~ f_mut320 in order
pred_fmut_sorted = [];

for i = 1:length(sorted_labels)
    label = sorted_labels{i};
    match_idx = [];
    for j = 1:length(dr)
        if ismember(label, dr(j).p)
            match_idx = j;
            break;
        end
    end

    if ~isempty(match_idx)
        pred_fmut_sorted = [pred_fmut_sorted; dr(match_idx).yExpSimu];
    else
        warning(['Label not found in dr: ', label]);
    end
end

wt_pred = dr(1).yExpSimu;  % WT prediction
pred = [wt_pred; pred_fmut_sorted];  % combine into 321 x 5

%% === Save all into final k2_res_new.mat ===
save('k2_res_new.mat', 'k2_res', 'pred', 'aa');


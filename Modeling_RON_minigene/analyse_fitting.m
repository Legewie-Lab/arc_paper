% Load data
load('hekrep1nowt.mat');
load('RON_Hek293_2_1000.mat');
load('C:\Users\Panajot\Nextcloud2\thesis\Writing\arc\arc_draft\corrections\submition\Codes\Coordinated-alternative-splicing-decisions-via-stepwise-exon-definition\RON data and fitting\Error_model_and_kd_fitting\datahekrep1.mat');

% Scale datar2 values to percentage
datar2(:, 2:6) = datar2(:, 2:6) .* 100;

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

% Assemble parameter tables
param_k1 = table(numk1, log10(k1_res .* log10(at.k1)));
param_k3 = table(numk3, log10(k3_res .* at.k3));

if moo == 1
    param_k2 = table(numk2, log10(k2_res .* log10(at.k2a)), log10(k2_res .* log10(at.k2b)));
else
    param_k2 = table(numk2, log10(k2_res .* log10(at.k2)));
end

% Split data and predictions
wt = data(1, :);
datak1 = data(2:lenk1 + 1, :);
datak2 = data(lenk1 + 2:lenk1 + lenk2 + 1, :);
datak3 = data(end - lenk3 + 1:end, :);
valk1 = pred(2:lenk1 + 1, :);
valk2 = pred(lenk1 + 2:lenk1 + lenk2 + 1, :);
valk3 = pred(end - lenk3 + 1:end, :);

% Fit curve for plotting
p = [-0.0894, 0.0914, 0.0140];
yw = linspace(0, 1);
vv = polyval(p, yw);
titre = ["Inclusion", "Skipping", "FULL IR", "1IR", "2IR"];

% Plot 1: Model vs Data with background fill
figure;
h1 = plot(data .* 100, pred .* 100, 'o');
hold on;

% Fill the background area (± model error)
fill([yw, fliplr(yw)] .* 100, [yw - vv, fliplr(yw + vv)] .* 100, ...
    [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'LineWidth', 0.1);

% Assign colors to each plot point
colorOrder = get(gca, 'ColorOrder');
for i = 1:min(size(data, 2), size(pred, 2))
    set(h1(i), 'Color', colorOrder(i, :));
end

% Plot the bisecting line
h2 = plot([0 1] .* 100, [0 1] .* 100, 'k-', 'LineWidth', 1.25);

% Create legend using valid graphics handles
legend([h1; h2], [titre, 'Bisecting line']);

xlabel('Data (%)');
ylabel('Model (%)');
xlim([0 100]);
ylim([0 100]);
title('Two-step model');
axis square;


% Plot 2: PSI vs Efficiency
figure;
da = [valk2(:, 1) ./ (valk2(:, 1) + valk2(:, 2)), 1 - sum(valk2(:, 3:5), 2)];
g = sortrows(da, 1) .* 100;

plot(datahekrep1{3}(:, 2) ./ (datahekrep1{3}(:, 2) + datahekrep1{3}(:, 3)) .* 100, ...
     (datahekrep1{3}(:, 2) + datahekrep1{3}(:, 3)) .* 100, 'go');
hold on;
plot(g(:, 1), g(:, 2), 'k', 'LineWidth', 1);
plot(datahekrep1{1}(:, 2) ./ (datahekrep1{1}(:, 2) + datahekrep1{1}(:, 3)) .* 100, ...
     (datahekrep1{1}(:, 2) + datahekrep1{1}(:, 3)) .* 100, 'k.', 'MarkerSize', 30);

xlabel('PSI (%)');
ylabel('Efficiency (%)');
xlim([0 100]);
ylim([40 100]);
title('Alternative exon mutation', 'FontSize', 18);
legend('Mutation data', 'Model fit', 'WT data', 'color', 'none');
axis square;












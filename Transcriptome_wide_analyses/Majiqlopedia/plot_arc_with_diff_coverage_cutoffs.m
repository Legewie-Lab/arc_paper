%%
clear
close all

coverage_cutoff = [10 20 50 80 100];

% coverage_cutoff = [10 20 40 80 160];

for i = 1:length(coverage_cutoff)
    arc{i} = load(sprintf('arc_coverage_cutoff_%d.mat',coverage_cutoff(i)));
end

% %%
% markers = {'o','s','*','^','+'};
% figure
% for i = 1:20
%     subplot(5,4,i)
%     for j = 1:length(coverage_cutoff)
%         plot(5:10:100,cellfun(@median,arc{j}.psi_eff_normal{i}),'b','Marker',markers{j});
%     hold on
%     end
%     for j = 1:length(coverage_cutoff)
%         plot(5:10:100,cellfun(@median,arc{j}.psi_eff_cancer{i}),'r','Marker',markers{j})
%     end
%     % plot(5:10:100,cellfun(@median,psi_eff_cancer{i}),'-o')
%     ylim([0 100])
%     title(arc{1}.inter_names{i},'FontWeight','normal')
% end

%
% psi_bins = 5:10:100;

% for j = 1:20
%     subplot()
    % for i = 1:length(psi_bins)
    %     len = length(arc{1}.psi_eff_normal{1}{i});
    %     cellfun(@(x)boxplot(x,num2cell(psi_bins(i)*ones(len,1))),arc{j}.psi_eff_normal{1}{i})
    %     hold on
    % end
% 
% end

% create long table
for i = 1:20

    for j = 1:5

        C_normal = arc{j}.psi_eff_normal{i};
        data_normal  = [];
        group_normal = [];
        
        C_cancer = arc{j}.psi_eff_cancer{i};
        data_cancer  = [];
        group_cancer = [];
        
        for k = 1:10
            
            data_normal  = [data_normal;  C_normal{k}];
            group_normal = [group_normal; repmat(k, numel(C_normal{k}), 1)];
            
            data_cancer  = [data_cancer;  C_cancer{k}];
            group_cancer = [group_cancer; repmat(k, numel(C_cancer{k}), 1)];
        
        end
        
        efficiency = data_normal;
        psi_bin = group_normal;
        % condition = repmat("normal",length(efficiency),1);
        condition = zeros(length(efficiency),1);
        tb_normal = table(efficiency,psi_bin,condition);
        
        
        efficiency = data_cancer;
        psi_bin = group_cancer;
        condition = ones(length(efficiency),1);
        tb_cancer = table(efficiency,psi_bin,condition);
        
        tb{j}{i} = [tb_normal;tb_cancer];
    end

end

% %%
% figure
% for i = 1: 10
%     for j = 1:5
%         subplot(10,5,(i-1)*5 + j)    
%         b=boxchart(tb{j}{i}.psi_bin,tb{j}{i}.efficiency,'GroupByColor',tb{j}{i}.condition);
%         b(1).MarkerSize = 1;
%         b(2).MarkerSize = 1;
%         b(1).MarkerColor = [200 220 255]/255; 
%         b(2).MarkerColor = [255, 204, 203]/255;
% 
%         xticks(1:10)
%         xticklabels('')
%         yticks(0:25:100)
%         ylim([-2 102])
%         set(gca,'TickDir','out')
% 
%         if i == 10
%             xticklabels({'0-10','10-20','20-30','30-40','40-50',...
%                 '50-60','60-70','70-80','80-90','90-100'})
%         end
% 
%         title(sprintf('\\color{blue}n_n=%d, \\color{red}n_c=%d', ...
%             sum(tb{j}{i}.condition == 0), ...
%             sum(tb{j}{i}.condition == 1)), ...
%                 'FontWeight','normal')
% 
%         if j == 1
%             ylabel(arc{1}.inter_names{i},'Rotation',0)
%         end
% 
%     end
% 
% end
% 
% width=21;
% height=width*1.4;
% 
% fig_position=[1,1,width,height];
% set(gcf,'units','centimeters');
% set(gcf,'paperunits','centimeters');
% set(gcf,'position',fig_position,'color','white');
% set(gcf,'paperposition',fig_position,'PaperSize',fig_position(3:4));
% 
% % save figure as .fig and .svg format
% fig_name = './results/arc_tissues_part1';  % change accordingly
% % fig_name = './results/arc_tissues_part1_exponential_cutoffs';  % change accordingly
% savefig(fig_name)
% print(fig_name,'-dsvg','-r1000','-vector')
% 
% %%
% figure
% for i = 1: 10
%     for j = 1:5
%         subplot(10,5,(i-1)*5 + j)    
%         b=boxchart(tb{j}{i+10}.psi_bin,tb{j}{i+10}.efficiency,'GroupByColor',tb{j}{i+10}.condition);
%         b(1).MarkerSize = 1;
%         b(2).MarkerSize = 1;
%         b(1).MarkerColor = [200 220 255]/255; 
%         b(2).MarkerColor = [255, 204, 203]/255;
% 
%         xticks(1:10)
%         xticklabels('')
%         yticks(0:25:100)
%         ylim([-2 102])
%         set(gca,'TickDir','out')
% 
%         if i == 10
%             xticklabels({'0-10','10-20','20-30','30-40','40-50',...
%                 '50-60','60-70','70-80','80-90','90-100'})
%         end
% 
%         % if j == 1
%         %     title(arc{1}.inter_names{i+10},'FontWeight','normal')
%         % end
% 
%         title(sprintf('\\color{blue}n_n=%d, \\color{red}n_c=%d', ...
%             sum(tb{j}{i+10}.condition == 0), ...
%             sum(tb{j}{i+10}.condition == 1)), ...
%                 'FontWeight','normal')
% 
%         if j == 1
%             ylabel(arc{1}.inter_names{i+10},'Rotation',0)
%         end
% 
% 
%     end
% 
% end
% 
% width=21;
% height=width*1.4;
% 
% fig_position=[1,1,width,height];
% set(gcf,'units','centimeters');
% set(gcf,'paperunits','centimeters');
% set(gcf,'position',fig_position,'color','white');
% set(gcf,'paperposition',fig_position,'PaperSize',fig_position(3:4));
% 
% % save figure as .fig and .svg format
% fig_name = './results/arc_tissues_part2';  % change accordingly
% % fig_name = './results/arc_tissues_part2_exponential_cutoffs';  % change accordingly
% savefig(fig_name)
% print(fig_name,'-dsvg','-r1000','-vector')

%%
% 

for i = 1:20
    tb_min = [];
    for j = 1:5
        tb_min_temp = tb{j}{i}(tb{j}{i}.psi_bin == 5 | tb{j}{i}.psi_bin == 6,:);
        cutoff = coverage_cutoff(j)*ones(size(tb_min_temp,1),1);
        % cutoff = (1.5*j-1)*ones(size(tb_min_temp,1),1);
        tb_min_temp = addvars(tb_min_temp,cutoff);
        tb_min = [tb_min;tb_min_temp];
    end
    tb_min.psi_bin = 0.5*ones(size(tb_min.psi_bin));
    tb_eff_min{i} = tb_min;

end

for i = 1:20
    tb_max = [];
    for j = 1:5
        tb_max_temp = tb{j}{i}(tb{j}{i}.psi_bin == 1 | tb{j}{i}.psi_bin == 10,:);
        cutoff = coverage_cutoff(j)*ones(size(tb_max_temp,1),1);
        % cutoff = (1.5*j-1)*ones(size(tb_max_temp,1),1);
        tb_max_temp = addvars(tb_max_temp,cutoff);
        tb_max = [tb_max;tb_max_temp];
    end
    tb_max.psi_bin = zeros(size(tb_max.psi_bin));
    tb_eff_max{i} = tb_max;
end

for i = 1:20
    tb_eff{i} = [tb_eff_min{i};tb_eff_max{i}];
end

%% bootstrapping median values
nboot = 1000;
for i = 1:20
    i
    for j = 1:5
        temp_normal = tb_eff{i}(tb_eff{i}.condition==0 & tb_eff{i}.cutoff==coverage_cutoff(j),:);
        temp_cancer = tb_eff{i}(tb_eff{i}.condition==1 & tb_eff{i}.cutoff==coverage_cutoff(j),:);
        log_median_ratio_normal{i}(:,j) = bootstrp(nboot,@calc_log_median_ratio,temp_normal);
        log_median_ratio_cancer{i}(:,j) = bootstrp(nboot,@calc_log_median_ratio,temp_cancer);
    end
end

function lmr = calc_log_median_ratio(df)

s = groupsummary(df,"psi_bin","median","efficiency");

if isscalar(unique(s.psi_bin))
    lmr = nan;
else
    % lmr = log10(s.median_efficiency(s.psi_bin == 0.5)/s.median_efficiency(s.psi_bin == 0));
    % lmr = s.median_efficiency(s.psi_bin == 0.5)/s.median_efficiency(s.psi_bin == 0);
    lmr = s.median_efficiency(s.psi_bin == 0.5) - s.median_efficiency(s.psi_bin == 0);

end

end


%
cutoff = repmat(coverage_cutoff,nboot,1);
cutoff = cutoff(:);
cutoff = [cutoff;cutoff];
condition = [zeros(length(cutoff)/2,1);ones(length(cutoff)/2,1)];

for i = 1:20

    log_median_ratio = [log_median_ratio_normal{i}(:); log_median_ratio_cancer{i}(:)];

    tb_lmr{i} = table(log_median_ratio,cutoff,condition);

end

%

cutoff_range = linspace(0, 110, 200)';

for i = 1:20

    ind_normal = tb_lmr{i}.condition == 0;
    ind_cancer = tb_lmr{i}.condition == 1;

    mdl_normal{i} = fitlm(tb_lmr{i}.cutoff(ind_normal),tb_lmr{i}.log_median_ratio(ind_normal));
    mdl_cancer{i} = fitlm(tb_lmr{i}.cutoff(ind_cancer),tb_lmr{i}.log_median_ratio(ind_cancer));

    % Get predictions with confidence intervals
    [eff_pred_normal{i}, eff_ci_normal{i}] = predict(mdl_normal{i}, cutoff_range, 'Alpha', 0.05, 'Simultaneous', false);
    [eff_pred_cancer{i}, eff_ci_cancer{i}] = predict(mdl_cancer{i}, cutoff_range, 'Alpha', 0.05, 'Simultaneous', false);

end

%%

for i = 1:20
    for j = 1:5
        j_normal = tb_lmr{i}.condition == 0 & tb_lmr{i}.cutoff == coverage_cutoff(j);
        j_cancer = tb_lmr{i}.condition == 1 & tb_lmr{i}.cutoff == coverage_cutoff(j);
        
        % [p_normal(i,j),h_normal(i,j)] = signtest(tb_lmr{i}.log_median_ratio(j_normal),0,'tail','left');
        % [p_cancer(i,j),h_cancer(i,j)] = signtest(tb_lmr{i}.log_median_ratio(j_cancer),0,'tail','left');
        [p_normal(i,j),h_normal(i,j)] = signrank(tb_lmr{i}.log_median_ratio(j_normal),0,'tail','left');
        [p_cancer(i,j),h_cancer(i,j)] = signrank(tb_lmr{i}.log_median_ratio(j_cancer),0,'tail','left');

    end
end

%%
figure
for i = 1:20
    subplot(4,5,i)
    b=boxchart(categorical(tb_lmr{i}.cutoff),tb_lmr{i}.log_median_ratio,'GroupByColor',tb_lmr{i}.condition);
    b(1).MarkerSize = 1;
    b(2).MarkerSize = 1;
    b(1).MarkerColor = [200 220 255]/255; 
    b(2).MarkerColor = [255, 204, 203]/255;

    % hold on
    % fill([cutoff_range/10; flipud(cutoff_range/10)], [eff_ci_normal{i}(:,1); flipud(eff_ci_normal{i}(:,2))], ...
    %  'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    % plot(cutoff_range/10, eff_pred_normal{i},'b')
    % 
    % fill([cutoff_range/10; flipud(cutoff_range/10)], [eff_ci_cancer{i}(:,1); flipud(eff_ci_cancer{i}(:,2))], ...
    %  'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    % plot(cutoff_range/10, eff_pred_cancer{i},'r')


    % xlim([0 110]/10)
    % xticks(coverage_cutoff/10)
    % xticks(unique(tb_eff_min{i}.cutoff))
    hold on
    yline(0)
    title(arc{1}.inter_names{i},'FontWeight','normal')
    xticklabels(split(num2str(coverage_cutoff)))
    ylim([-91 20])
    set(gca,'TickDir','out')

    if i == 18
        xlabel('Coverage cutoff')
    end

    if i == 6
        ylabel('\Deltaefficiency between intermediate and extreme PSIs (%)')
    end

    if sum(h_normal(i,:)) < 5
        ind_null = find(h_normal(i,:) == 0);
        for ii = 1:length(ind_null)
            text(categorical(coverage_cutoff(ind_null(ii))),15,'ns')
        end
    end

    if sum(h_cancer(i,:)) < 5
        ind_null = find(h_cancer(i,:) == 0);
        for ii = 1:length(ind_null)
            text(categorical(coverage_cutoff(ind_null(ii))),15,'ns')
        end
    end

end

width=21;
height=width/1.2;

fig_position=[1,1,width,height];
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'position',fig_position,'color','white');
set(gcf,'paperposition',fig_position,'PaperSize',fig_position(3:4));

% save figure as .fig and .svg format
fig_name = './results/delta_eff_at_mid_psi';  % change accordingly
savefig(fig_name)
print(fig_name,'-dsvg','-r1000','-vector')

%% fit median values directly without bootstrapping
% %%
% 
% for i = 1:20
%     tb_summary_eff{i} = groupsummary(tb_eff{i},["condition","cutoff","psi_bin"],"median","efficiency");
% end
% 
% 
% %%
% for i = 1:20
%     tb_summary_eff_min{i} = groupsummary(tb_eff_min{i},["condition","cutoff"],"median","efficiency");
%     tb_summary_eff_max{i} = groupsummary(tb_eff_max{i},["condition","cutoff"],"median","efficiency");
% 
%     log_median_ratio_normal{i} = log10(tb_summary_eff_min{i}.median_efficiency(tb_summary_eff_min{i}.condition == 0) ...
%                                      ./tb_summary_eff_max{i}.median_efficiency(tb_summary_eff_max{i}.condition == 0));
%     log_median_ratio_cancer{i} = log10(tb_summary_eff_min{i}.median_efficiency(tb_summary_eff_min{i}.condition == 1) ...
%                                      ./tb_summary_eff_max{i}.median_efficiency(tb_summary_eff_max{i}.condition == 1));
% 
%     for j = 1:5
%         p_normal{i}(j,1) = ranksum(tb_eff_max{i}.efficiency(tb_eff_max{i}.cutoff == coverage_cutoff(j) ...
%                                  & tb_eff_max{i}.condition == 0), ...
%                                    tb_eff_min{i}.efficiency(tb_eff_min{i}.cutoff == coverage_cutoff(j) ...
%                                  & tb_eff_min{i}.condition == 0));
%         p_cancer{i}(j,1) = ranksum(tb_eff_max{i}.efficiency(tb_eff_max{i}.cutoff == coverage_cutoff(j) ...
%                                  & tb_eff_max{i}.condition == 1), ...
%                                    tb_eff_min{i}.efficiency(tb_eff_min{i}.cutoff == coverage_cutoff(j) ...
%                                  & tb_eff_min{i}.condition == 1));
% 
%     end
% end
% 
% for i = 1:20
%     subplot(4,5,i)
%     scatter(coverage_cutoff,log_median_ratio_normal{i})
%     hold on
%     scatter(coverage_cutoff,log_median_ratio_cancer{i})
%     ylim([-0.15 0])
% end
    

% cutoff_range = linspace(0, 110, 200)';
% 
% for i = 1:20
% 
%     tb_summary = groupsummary(tb_eff_min{i},["condition","cutoff"],"median","efficiency");
% 
%     mdl_normal{i} = fitlm(tb_summary.cutoff(tb_summary.condition == 0),tb_summary.median_efficiency(tb_summary.condition == 0));
%     mdl_cancer{i} = fitlm(tb_summary.cutoff(tb_summary.condition == 1),tb_summary.median_efficiency(tb_summary.condition == 1));
% 
%     % Get predictions with confidence intervals
%     [eff_pred_normal{i}, eff_ci_normal{i}] = predict(mdl_normal{i}, cutoff_range, 'Alpha', 0.05, 'Simultaneous', false);
%     [eff_pred_cancer{i}, eff_ci_cancer{i}] = predict(mdl_cancer{i}, cutoff_range, 'Alpha', 0.05, 'Simultaneous', false);
% end
% 
% for i = 1:20
%     k_normal(i)=mdl_normal{i}.Coefficients.Estimate(end);
%     k_cancer(i)=mdl_cancer{i}.Coefficients.Estimate(end);
% end
% %
% figure
% bar([k_normal;k_cancer]')
% 
% %
% figure
% for i = 1:20
%     subplot(5,4,i)
%     b=boxchart(tb_eff_min{i}.cutoff/10,tb_eff_min{i}.efficiency,'GroupByColor',tb_eff_min{i}.condition);
%     b(1).MarkerSize = 1;
%     b(2).MarkerSize = 1;
%     b(1).MarkerColor = [200 220 255]/255; 
%     b(2).MarkerColor = [255, 204, 203]/255;
% 
%     hold on
%     fill([cutoff_range/10; flipud(cutoff_range/10)], [eff_ci_normal{i}(:,1); flipud(eff_ci_normal{i}(:,2))], ...
%      'b', 'EdgeColor', 'none', 'FaceAlpha', 0.25)
%     plot(cutoff_range/10, eff_pred_normal{i},'b')
% 
%     fill([cutoff_range/10; flipud(cutoff_range/10)], [eff_ci_cancer{i}(:,1); flipud(eff_ci_cancer{i}(:,2))], ...
%      'r', 'EdgeColor', 'none', 'FaceAlpha', 0.25)
%     plot(cutoff_range/10, eff_pred_cancer{i},'r')
% 
% 
%     xlim([0 110]/10)
%     xticks(coverage_cutoff/10)
%     % xticks(unique(tb_eff_min{i}.cutoff))
%     xticklabels(split(num2str(coverage_cutoff)))
%     ylim([-5 105])
%     set(gca,'TickDir','out')
% 
% end


%%
% %
% figure
% for i = 1:20
%     subplot(5,4,i)
%     for j = 1:length(coverage_cutoff)
%         plot(5:10:100,cellfun(@median,arc{j}.psi_eff_cancer{i}),'-o')
%     hold on
%     end
%     % plot(5:10:100,cellfun(@median,psi_eff_cancer{i}),'-o')
%     ylim([0 100])
% end


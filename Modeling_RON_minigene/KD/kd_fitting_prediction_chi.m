%% kd_fitting_prediction_Jay.m
% Fit and summarize RBP knockdown predictions using the k2 perturbation model.
% This version reports epsilon parameters and exports WT_epsilon_parameters.xlsx.

%% Loading data
clear all
close all


load('datahekrep1nowt.mat')
load('datamcf7avgnowt.mat')
load('datamcf7kdavgnowt.mat')
load('dataprpf6nowt.mat')
load('datasrsf2nowt.mat')
load('datasmu1nowt.mat')
load('datakontrollenowt.mat')
load('datapuf60nowt.mat')
load('datahnrnph1nowt.mat')

load('dataprpf6.mat')
load('datasrsf2.mat')
load('datasmu1.mat')
load('datakontrolle.mat')
load('datapuf60.mat')
load('datahnrnph1.mat')



% error model
% load('error_fitpred.mat')

load('error_fitpred_new.mat')
error_fitpred = error_fitpred_new;

% perturbation factor
load('k2_res.mat')
load('aa.mat')


% load('k2_res_new.mat')
% load('k2_res_new_from_best_control.mat')
% load('k2_res_new_from_best_control_SE.mat')

load('k2_res_new_from_best_control_SE0.mat')
% if istable(aa)
%     aa = table2array(aa);
% end

%% fitting parameters and data matrix ---
errorwt=error_fitpred(1,:);

num=3;

res=repmat(aa,size(k2_res));
res(:,2:3)=aa(2:3).*k2_res;
% res=[aa;res];

% WT + Exon2
dkontrolle=[datakontrollenowt{1};datakontrollenowt{3}];

% data matrix
% WT + Exon2
alldat=[{[dataprpf6nowt{1};dataprpf6nowt{3}]};{[datapuf60nowt{1};datapuf60nowt{3}]};{[datahnrnph1nowt{1};datahnrnph1nowt{3}]};{[datasmu1nowt{1};datasmu1nowt{3}]};{[datasrsf2nowt{1};datasrsf2nowt{3}]};{[datamcf7avgnowt{1};datamcf7avgnowt{3}]}];
% WT
alldatwt=[{[dataprpf6nowt{1}]};{[datapuf60nowt{1}]};{[datahnrnph1nowt{1}]};{[datasmu1nowt{1}]};{[datasrsf2nowt{1}]};{[datamcf7avgnowt{1}]}];


% intersection index of mutant vs control.
cik2=[];
cibk2=[];

cik2wt=[];
cibk2wt=[];

for i=1:length(alldat)
    [Ca,iak2,ibk2]=intersect(alldat{i}(:,1),dkontrolle(:,1));
    cik2=[cik2;{iak2}];
    cibk2=[cibk2;{ibk2}];

    [Cawt,iak2wt,ibk2wt]=intersect(alldatwt{i}(:,1),dkontrolle(:,1));
    cik2wt=[cik2wt;{iak2wt}];
    cibk2wt=[cibk2wt;{ibk2wt}];
end

indx=[cik2,cibk2];
indxwt=[cik2wt,cibk2wt];

%% Multi-model fitting (only WT)
tic
sol=multik2lsomikronfindercommit(alldatwt,indxwt,aa,error_fitpred(1,:));
toc

% --- Extract WT epsilon (b) parameters 
RBP_names = {'PRPF6','PUF60','HNRNPH1','SMU1','SRSF2'};
num_RBP = 5;        
num_cond = 15;

B_table = zeros(num_cond, num_RBP);

for j = 1:num_RBP
    B_table(:,j) = sol{j,2};   % extract b values
end

% column
indication = {'k1','k2a','k2b','k2','k3',...
              'k1&k2a','k1&k2b','k1&k2','k1&k3',...
              'k2a&k3','k2b&k3','k2&k3',...
              'k1&k2a&k3','k1&k2b&k3','kon'};

% turn to table
T = table(indication', B_table(:,1), B_table(:,2), ...
          B_table(:,3), B_table(:,4), B_table(:,5), ...
          'VariableNames', {'Indication','PRPF6','PUF60','HNRNPH1','SMU1','SRSF2'});

disp(T)

writetable(T, 'WT_epsilon_parameters.xlsx')















%% Predict the isoform frequency of mutants using the best-fit parameter

% sol=multik2lsomikronfindercommit(alldatwt,indx,res,num,errorek2);
tic
% val=[];
ksol=[];
for j=1:length(sol)
    val=[];

    for cond=1:length(sol{1})
        y=multiomikroncalculatorcalccommit(sol{j,2}(cond),res(indx{j,2},:),alldat{j}(indx{j,1},2:6),0,cond,1);

        val=[val;{y}]; 
%         val(cond,:)={y};
%         cond
    end
    ksol=[ksol;{val}];
end
toc

%%  chi-square fit score
for j=1:length(ksol)
    chi=[];
    for cond=1:15
        divf=sum((((alldat{j}(indx{j,1}(2:end,:),2:6)-ksol{j}{cond}(2:end,:))./error_fitpred(indx{j,2}(2:end,:),2:6)).^2)./(chi2inv(0.95,numel(alldat{j}(indx{j,1}(2:end,:),2:6)))),'all');
        chi=[chi;{divf}];

    end
    ksol{j,2}=chi;
end

% Model/Condition Naming
legName=[{'PRPF6'},{'PUF60'},{'HNRNPH1'},{'SMU1'},{'SRSF2'},{'MCF7'}];
    % indication=['spliceosome binding perturbation in exon1 ',"U2' of exon 2 binding perturbation ","U1' of exon 2 binding perturbation  ",'spliceosome binding perturbation in exon2','spliceosome binding perturbation in exon3',"spliceosome binding perturbations in exon1 and U2' of exon 2",'spliceosome binding perturbations in exon1 and exon2','spliceosome binding perturbation in exon1 and exon3',"spliceosome binding perturbation in exon3 and U1' of exon 2",'spliceosome binding perturbation in exon2 and exon3','spliceosome binding perturbation in all three exons'];
    % xvalues = {'MCF7 k2','PRPF6 k2a','SRSF2 kon','SMU1 k2&k3','PUF60 k2b','HRNRPH1 k2'};
indication={'k1 ',"k2a ","k2b","k2  ",'k3','k1&k2a',"k1&k2b" ,"k1&k2",'k1&k3','k2a&k3','k2b&k3',"k2&k3",'k1&k2a&k3','k1&k2b&k3','kon'};


% Fitting quality results
solsd=[];
soleps=[];
solmat=[];
for i=1:size(sol,1);
   solsd(:,i)=sol{i,1};
   soleps(:,i)=sol{i,2};

   solmat(:,i)=cellfun(@double, ksol{i,2});
end




%% Heatmap
xvalue_wts = {'PRPF6  control','PUF60  control','HRNRPH1  control','SMU1  control','SRSF2  control','MCF7 control'};
xvalue_nor = {'PRPF6 prediction','PUF60 prediction','HRNRPH1 prediction','SMU1 prediction','SRSF6 prediction','MCF7 prediction'};
xvalue_na = {'PRPF6','PUF60','HRNRPH1','SMU1','SRSF2'};

newXvalues=[];

for i=1:5
    newXvalues=[newXvalues,{xvalue_wts{i}},{xvalue_nor{i}}];
end
titre=["Inclusion" "Skipping" "FullIR" "FirstIR" "SecondIR"];

combinedMatrix = [solsd, solmat];

% Assuming A and B are already defined as 15x6 matrices

% Determine the size of the matrices A and B
[rows, cols] = size(solsd);

% Preallocate matrix C
C = zeros(rows, 2*cols);

% Populate matrix C with alternating values from A and B
for col = 1:cols
    % Compute the column indices in C for A and B
    colIndices = (2*col-1):(2*col);
    
    % Assign values from A and B to the corresponding columns in C
    C(:, colIndices) = [solsd(:, col), solmat(:, col)];
end

% %%%%%%%%%%%%%%%%%  Heatmap 1: Prediction error (white)   %%%%%%%%%%%%%%%%% 
figure
h= heatmap(xvalue_na,indication,(round(soleps(:,1:5),2)));

h.Colormap = [1 1 1];  % White color

% %%%%%%%%%%%%%%%%%  Heatmap 2: Fitted vs. Predicted chi² values ​​(color heatmap)   %%%%%%%%%%%%%%%%% 
figure
h= heatmap(newXvalues,indication,(round(C(:,1:10),2)));
h.ColorScaling = 'scaledcolumns';
h.Colormap = hot;

% %%%%%%%%%%%%%%%%% Figure 5B    %%%%%%%%%%%%%%%%% 
% log10 transform
log_chi2_fit  = log10(solsd);     
log_chi2_pred = log10(solmat);   

% interleave fitting and prediction columns
[rows, cols] = size(log_chi2_fit);  
C = zeros(rows, 2 * cols);

for col = 1:cols
    C(:, 2*col - 1) = log_chi2_fit(:, col);    % fitting
    C(:, 2*col)     = log_chi2_pred(:, col);   % prediction
end

% labels
xlabels = {'PRPF6 fitting','PRPF6 prediction','PUF60 fitting','PUF60 prediction',...
           'HNRNPH1 fitting','HNRNPH1 prediction','SMU1 fitting','SMU1 prediction',...
           'SRSF2 fitting','SRSF2 prediction','MCF7 fitting','MCF7 prediction'};

% plot heatmap
figure;
h = heatmap(xlabels, indication, C);

% Colormap: grayscale from -2 (good) to +1 (bad)
h.Colormap = bone;
h.ColorLimits = [-1.5 2.5];
h.CellLabelFormat = '%.2f';
% h.CellLabelColor = 'none';   
title('Figure 5B-style: log_{10}(\chi^2 / \chi^2_{crit})');


figure5B_data.solsd       = solsd;
figure5B_data.solmat      = solmat;
figure5B_data.C           = C;
figure5B_data.indication  = indication;
figure5B_data.xlabels     = xlabels;
figure5B_data.metric      = 'RMSE';
figure5B_data.transform   = 'log10';
figure5B_data.date        = datestr(now);

save('Figure5B_data.mat', 'figure5B_data');


% -- Plot 2
S = load('Figure5B_data.mat');
data = S.figure5B_data;

figure;
h = heatmap(data.xlabels, data.indication, data.C);

h.Colormap = bone;
h.CellLabelFormat = '%.2f';

title('Chi2 comparison: Fitting vs Prediction');


%% Data vs. Fitted Plot (Panel in Figure 5C)
% Find the optimal model (minimum chi²) for each RBP.
lochi=[];
for i=1:length(sol)
    [mini,idex]=min(sol{i,1});
    lochi=[lochi,idex];
end

%%%%%%%%%%%%%%%%% Data vs. Fitted %%%%%%%%%%%%%%%%%
figure
for j=1:6

    nexttile
    plot(alldat{j}(indx{j,1},2:6),ksol{j}{lochi(j)},'o'),hold on
    plot([0 1],[0,1],'k')
    % set(gcf,'Color','k'),set(gca,'XColor','w'),set(gca,'YColor','w','ZColor','w','Color','k')

    xlabel('data'),ylabel('prediction'),title([legName{j},indication(lochi(j)),num2str(sol{j,1}(lochi(j)))],'Color','k')

%     sgtitle(legName{j})
end
% legend(titre,'Color','k','TextColor','w')

%%%%%%%%%%%%%%%%%  PSI vs Efficiency %%%%%%%%%%%%%%%%%
figure
plot(dkontrolle(:,2)./(dkontrolle(:,2)+dkontrolle(:,3)),1-sum(dkontrolle(:,4:6),2),'o',pred(:,1)./(pred(:,1)+pred(:,2)),1-sum(pred(:,3:5),2),'.')
hold on
plot(dkontrolle(1,2)./(dkontrolle(1,2)+dkontrolle(1,3)),1-sum(dkontrolle(1,4:6),2),'g*',pred(1,1)./(pred(1,1)+pred(1,2)),1-sum(pred(1,3:5),2),'k+')
xlabel('PSI'),ylabel('Efficiency'),title('datavkfitting kontrolle exon2')
legend('data','fitting','wt data','wtfitting')

% p=[-0.2316,0.2420,0.0135];
de=linspace(0,1);
p=[-0.1323    0.1354    0.0137];%polifit mean std control data
sd=polyval(p,de);
% sp=polyval(p,pred);
% plot(de,se,'.')

%%%%%%%%%%%%%%%%%  Five isoform fitted values ​​vs. data point scatter plots    %%%%%%%%%%%%%%%%%
figure
plot(dkontrolle(:,2:6).*100,pred.*100,'o',de.*100,(de-sd).*100,'k',de.*100,(de+sd).*100,'k',[0 1].*100, [0 1].*100,'k')
xlabel('data'),ylabel('fitting'),xlim([0 1].*100),ylim([0 1].*100),legend(titre),title('datavkfitting kontrolle exon2')

titre=["Inclusion" "Skipping" "FullIR" "FirstIR" "SecondIR"];
% p=[-0.2316 0.2420 0.0135];
% shaded area
yw=linspace(0,1);
vv=polyval(p,yw,3);

correr=[];


%%%%%%%%%%%%%%%%%  Figure S5 (middel): scatter plot of Model predictions vs. experimental values   %%%%%%%%%%%%%%%%%
figure
for j=1:5
subplot(3,5,j+5)
    % nexttile
    plot(alldat{j}(indx{j,1},2:6) .*100,ksol{j}{lochi(j)} .*100,'o'),hold on
    plot( yw.*100, yw.*100 + vv.*100, 'k-', yw.*100, yw.*100 - vv.*100, 'k-',[0 100], [0, 100], 'k--', 'LineWidth', 1.5)
    fill([yw.*100, fliplr(yw.*100)], [yw.*100 - vv.*100, fliplr(yw.*100 + vv.*100)], [0.7 0.7 0.7] ,'FaceAlpha', 0.3)

    pointsInside = inpolygon(reshape((alldat{j}(indx{j, 1}, 2:6)).', 1, []),reshape((ksol{j}{lochi(j)}).', 1, []), [yw.*100, fliplr(yw.*100)], [yw.*100 - vv.*100, fliplr(yw.*100 + vv.*100)]);

    xlabel('data'),ylabel('Model fitting'),xlim([0 100]),ylim([0 100])
    numPointsInside = (sum(pointsInside,'all')./length(alldat{1}(:))).*100;
    oo=rcorr(reshape(alldat{j}(indx{j, 1}, 2:6).', 1, []),reshape(ksol{j}{lochi(j)}.', 1, []));
    correr=[correr;oo];
    % text(0.1,0.8,['points insidel=',num2str(round(numPointsInside,1)),'%'])
    text(0.1.*100,0.7.*100,['R=',num2str((oo).*100),'%'])

    xlabel('Data'),ylabel('Prediction'),title([legName{j}])
    axis square


%     sgtitle(legName{j})
end
legend([titre])
% legend(titre,'Color','k','TextColor','w'
% print('C:\Users\Panajot\Nextcloud2\thesis\Writing\arc\arc_figures\supplement\figure_smult', '-dsvg', '-r2000', '-painters')


%%%%%%%%%%%%%%%%%  Figure 5C: PSI vs Efficiency    %%%%%%%%%%%%%%%%% 
figure
for j=[2,4]

    nexttile
    plot(alldat{j}(indx{j,1},2:6).*100,ksol{j}{lochi(j)}.*100,'o'),hold on
    plot([0 1].*100, [0, 1].*100, 'k')
    fill([yw, fliplr(yw)].*100, [yw - vv, fliplr(yw + vv)].*100, [0.7 0.7 0.7], 'FaceAlpha', 0.4,'LineWidth', 0.1)

    pointsInside = inpolygon(reshape(alldat{j}(indx{j, 1}, 2:6).', 1, []),reshape(ksol{j}{lochi(j)}.', 1, []), [yw, fliplr(yw)], [yw - vv, fliplr(yw + vv)]);

    xlabel('Data(%)'),ylabel('Model(%)'),xlim([0 1].*100),ylim([0 1].*100)
    numPointsInside = (sum(pointsInside,'all')./length(alldat{1}(:))).*100;
    oo=rcorr(reshape(alldat{j}(indx{j, 1}, 2:6).', 1, []),reshape(ksol{j}{lochi(j)}.', 1, []));
    correr=[correr;oo];
%     text(0.1,0.8,['points inside=',num2str(round(numPointsInside,1)),'%'],'FontSize',13)
%     text(0.1,0.7,['R=',num2str((oo).*100),'%'],'FontSize',13)

    xlabel('Data(%)'),ylabel('Model(%)'),title(legName{j})
    axis square

end
legend([titre,'Bisecting line'],'color','none')
% print('C:\Users\Panajot\Nextcloud2\thesis\Writing\arc\arc_figures\figure_4c', '-dsvg', '-r2000', '-painters')


%% data & Plot
alldat=[{[dataprpf6{1};dataprpf6{3}]};{[datapuf60{1};datapuf60{3}]};{[datahnrnph1{1};datahnrnph1{3}]};{[datasmu1{1};datasmu1{3}]};{[datasrsf2{1};datasrsf2{3}]};{[datamcf7avgnowt{1};datamcf7avgnowt{3}]}];

data_co=[{vertcat(dataprpf6{:})},{vertcat(datapuf60{:})},{vertcat(datahnrnph1{:})},{vertcat(datasmu1{:})},{vertcat(datasrsf2{:})}];

contro_matx=vertcat(datakontrolle{:});
%%%%%%%%%%%%%%%%%  Figure S5: PSI vs Efficiency  (Top + btm) %%%%%%%%%%%%%%%%% 
figure
best_fi={};
best_pi={};

for j=1:5
        % nexttile
        subplot(2,5,j+5)
        plot(alldat{j}(:,2)./(alldat{j}(:,2)+alldat{j}(:,3)),1-sum(alldat{j}(:,4:6),2),'o',ksol{j}{lochi(j)}(:,1)./(ksol{j}{lochi(j)}(:,1)+ksol{j}{lochi(j)}(:,2)),1-sum(ksol{j}{lochi(j)}(:,3:5),2),'.'),hold on
        plot(alldat{j}(1,2)./(alldat{j}(1,2)+alldat{j}(1,3)),1-sum(alldat{j}(1,4:6),2),'k*',ksol{j}{lochi(j)}(1,1)./(ksol{j}{lochi(j)}(1,1)+ksol{j}{lochi(j)}(1,2)),1-sum(ksol{j}{lochi(j)}(1,3:5),2),'c+')
        best_fi{j}=[ksol{j}{lochi(j)}];
        best_pi{j}=[ksol{j}{lochi(j)}(:,1)./(ksol{j}{lochi(j)}(:,1)+ksol{j}{lochi(j)}(:,2)),1-sum(ksol{j}{lochi(j)}(:,3:5),2)];
        xlabel('PSI'),ylabel('Efficiency'),title([legName{j},indication(lochi(j)),num2str(sol{j,1}(lochi(j)))])
        axis square
end
legend('data','fitting','wt data','wt fitting')


best_sort_pi = cellfun(@(x) sortrows(x,1), best_pi, 'UniformOutput', false);


% figure
for j=1:5
        % nexttile
        subplot(3,5,j+10)
        plot((alldat{j}(:,2)./(alldat{j}(:,2)+alldat{j}(:,3))).*100,(1-sum(alldat{j}(:,4:6),2)).*100,'o',best_sort_pi{j}(:,1).*100,best_sort_pi{j}(:,2).*100),hold on
       
        plot((alldat{j}(1,2)./(alldat{j}(1,2)+alldat{j}(1,3))).*100,(1-sum(alldat{j}(1,4:6),2)).*100,'k*',(ksol{j}{lochi(j)}(1,1)./(ksol{j}{lochi(j)}(1,1)+ksol{j}{lochi(j)}(1,2))).*100,(1-sum(ksol{j}{lochi(j)}(1,3:5),2)).*100,'c+')
        xlabel('PSI (%)'),ylabel('Efficiency (%)'),title(legName{j})
        axis square
end
legend('Data','Prediction','Wt Data','Wt Fitting')

% print('C:\Users\Panajot\Nextcloud2\thesis\Writing\arc\arc_figures\supplement\figure_s8', '-dsvg')


ihek=[];
ikd=[];

for i=1:length(data_co);
    [Cd,ak,bk]=intersect(contro_matx(:,1),data_co{i}(:,1));
    ihek=[ihek;{ak}];
    ikd=[ikd;{bk}];
end

% figure
% tiledlayout(1,5,'TileSpacing',"tight",'Padding','tight')

for i=1:length(data_co)
    subplot(3,5,i)
    plot(contro_matx(ihek{i},2:6).*100,data_co{i}(ikd{i},2:6).*100,'o',[0,1].*100,[0,1].*100,'k')
    xlabel('Control (%)'),ylabel('KD (%)'),title(legName{i}),xlim([0 100]),ylim([0 100])
    axis square
end
legend(titre)

%%%%%%%%%%%%%%%%%  Figure 5D: PSI vs. Efficiency  %%%%%%%%%%%%%%%%%
figure
for j=[2,4]
        nexttile
        plot([alldat{j}(:,2)./(alldat{j}(:,2)+alldat{j}(:,3))].*100,[1-sum(alldat{j}(:,4:6),2)].*100,'bo'),hold on
        da=[ksol{j}{lochi(j)}(:,1)./(ksol{j}{lochi(j)}(:,1)+ksol{j}{lochi(j)}(:,2)),1-sum(ksol{j}{lochi(j)}(:,3:5),2)];
        g=sortrows(da, 1).*100;
        plot(g(:,1),g(:,2),'r'),hold on
        plot([alldat{j}(1,2)./(alldat{j}(1,2)+alldat{j}(1,3))].*100,[1-sum(alldat{j}(1,4:6),2)].*100,'k.','MarkerSize', 30),hold on
        plot([ksol{j}{lochi(j)}(1,1)./(ksol{j}{lochi(j)}(1,1)+ksol{j}{lochi(j)}(1,2))].*100,[1-sum(ksol{j}{lochi(j)}(1,3:5),2)].*100,'r.','MarkerSize',30)
%         xlabel('PSI(%)'),ylabel('Efficiency(%)'),title([legName{j},indication(lochi(j)),num2str(sol{j,1}(lochi(j)))])
        xlabel('PSI(%)'),ylabel('Efficiency(%)'),title(legName{j})

        axis square
%         set(gca, 'FontSize', 14);  % Adjust the font size as needed
%         set(gca, 'LineWidth', 1.5); 


end
legend('data','prediction','wt (data)','wt (fit)','color','none');



%% Check
j = 2;   % PUF60 (PRPF6=1, PUF60=2, HNRNPH1=3, ...)

new_chi2 = zeros(15,1);

% --- KD experimental data ---
data = alldat{j}(indx{j,1}, 2:6);    % N × 5

% --- sigma (must match fitting sigma) ---
sigma = error_fitpred(indx{j,2}, 2:6);   % N × 5

for cond = 1:15

    % --- prediction using fitted parameter b ---
    pred = multiomikroncalculatorcalccommit( ...
        sol{j,2}(cond), ...           % fitted b
        res(indx{j,2},:), ...         % kinetic parameters
        data, ...
        0, cond, sigma);

    % --- chi-square ---
    r = (pred - data) ./ sigma;
    new_chi2(cond) = sum(r(:).^2);

end

df = numel(data)-1;  % or numel(data)-1
chi2crit = chi2inv(0.95, df);

new_chi2_norm = new_chi2 ./ chi2crit;
log10(new_chi2_norm)


%% Function
% multik2lsomikronfindercommit(wtkd,inde,aa,error);
   % fitting using WT(15 model)
   % The structure of sol : each RBP has 15 b and a corresponding normalized chi².

% multiomikroncalculatorcalccommit(b,aa,output,fit,cond,error)
   % Iterate through each RBP and predict KD data using the best-fit parameters for its 15 cond models.
   % obtain ksol{j}{cond}, which represents the isoform distribution predicted by the model under each mechanism assumption.

% mech_reverse_syndraftallplotcommit(t,x,k1,k2a,k2b,k3,kret,k4,k5a,k5b,k6,s,kspli,kincl,kskip,kdr1,kdr2,num)
   % The actual dynamic differential equations (22 variables)
   % Simulate the dynamic process of splicing from the initial state to the steady state.


function sol = multik2lsomikronfindercommit(wtkd,inde,aa,error);


    options=optimset('display','off');
    
    
    sol=[];
    z=[];
    ls=[];
    for j=1:length(wtkd)
        z=[];
        ls=[];
    %     tol=[];
        parfor cond =1:15;
            
            [w,pp]=lsqnonlin(@(b) multiomikroncalculatorcalccommit(b,aa(inde{j,2},:),wtkd{j}(inde{j,1},2:6),1,cond,error(inde{j,2},2:6)),1,0,20,options);
            z=[z;w];
            ls=[ls;pp./(chi2inv(0.95,numel(wtkd{j})))];
            
    
    
        end
    
        sol=[sol;{ls},{z}];
        
    end
end




function [y,o] = multiomikroncalculatorcalccommit(b,aa,output,fit,cond,error)
    
    % Y2=[];
    Y1=[];
    
    
    k1=aa(:,1);
    k2a=aa(:,2);
    k2b=aa(:,3);
    k3=aa(:,4);
    kret=aa(:,5);
    % kd:nt=aa(5);
    k4=aa(:,6);
    k5a=aa(:,7);
    k5b=aa(:,8);
    k6=aa(:,9);
    kspli=aa(:,10);
    kincl=aa(:,11);
    kskip=aa(:,12);
    kdr1=aa(:,13);
    kdr2=aa(:,14);
    s=aa(:,15);
    cop=1;
    
    
    if cond==1
        k1=k1.*b;
    elseif cond==2;
        k2a=k2a.*b;
    elseif cond==3;
        k2b=k2b.*b;
    elseif cond==4;
        k2a=k2a.*(b./1);
        k2b=k2b.*(b./1);
    elseif cond==5;
        k3=k3.*b;
    elseif cond==6;
        k1=k1.*b;
        k2a=k2a.*b;
    elseif cond==7;
        k1=k1.*b;
        k2b=k2b.*(b./1);
    elseif cond==8;
        k1=k1.*b;
        k2a=k2a.*(b./1);
        k2b=k2b.*(b./1);
    elseif cond==9;
        k1=k1.*b;
        k3=k3.*b;
    elseif cond==10;
        k2a=k2a.*b;
        k3=k3.*b;
    elseif cond==11;
        k2b=k2b.*b;
        k3=k3.*b;
    elseif cond==12;
        k2a=k2a.*(b./1);
        k2b=k2b.*(b./1);
        k3=k3.*b;
    elseif cond==13;
        k1=k1.*b;
        k2a=k2a.*(b./1);
        k3=k3.*b;
    elseif cond==14;
        k1=k1.*b;
        k2b=k2b.*(b./1);
        k3=k3.*b;
    elseif cond==15;
        k1=k1.*b;
        k2a=k2a.*(b./1);
        k2b=k2b.*(b./1);
        k3=k3.*b;
    
    end
        
        
        % for cond=1:size(output,1)
        
        
        ve=22;
        x0= zeros(1,ve);
        % x0(1)=1;
        options = odeset('NonNegative',1:ve,'RelTol',1e-8,'AbsTol',[1e-8]);
        
        t=linspace(0,1000000000, 100);
        %     f=@(t,x) mech_reverse_ret1(t,x0,b,k2(i));
    for i=1:size(output,1)
        [T,Y]=ode15s(@(t,x) mech_reverse_syndraftallplotcommit(t,x,k1(i,:),k2a(i,:),k2b(i,:),k3(i,:),kret(i,:),k4(i,:),k5a(i,:),k5b(i,:),k6(i,:),s(i,:),kspli(i,:),kincl(i,:),kskip(i,:),kdr1(i,:),kdr2(i,:)),t,x0,options);
        Y1=[Y1 ;Y(end,:)];
    end
    
    incl=Y1(:,13);
    skip=Y1(:,14);
    fullIR=[Y1(:,15)+sum(Y1(:,1:12),2)];
    fIR=Y1(:,16)+Y1(:,17)+Y1(:,20);
    seIR=Y1(:,18)+Y1(:,19)+Y1(:,21);
    su=[incl+skip+fullIR+fIR+seIR];
    Y1=[incl./su skip./su fullIR./su fIR./su seIR./su];
    
    
    if fit==1
        y = ((Y1-output)./error);
    
    else
        y = Y1 ; 
        o=[k1,k2a,k2b,k3,kret,k4,k5a,k5b,k6,kspli,kincl,kskip,kdr1,kdr2,s];
    end
end


function dxdt = mech_reverse_syndraftallplotcommit(t,x,k1,k2a,k2b,k3,kret,k4,k5a,k5b,k6,s,kspli,kincl,kskip,kdr1,kdr2,num)

    P000=x(1);
    P100=x(2);
    P0a0=x(3);
    P010=x(4);
    P001=x(5);
    P1a0=x(6); 
    P110=x(7); 
    P101=x(8); 
    P0a1=x(9);
    P011=x(10);
    P1a1=x(11);
    P111=x(12);
    incl=x(13);
    skip=x(14);
    ret=x(15);
    P011x=x(16);
    P111x=x(17);
    Px110=x(18);
    Px111=x(19);       
    ret1=x(20);
    ret2=x(21);
    deg=x(22);
    


	dP000dt = s+k4*P100+k5a*P0a0+k6*P001-(k1+k2a+k3+kret)*P000;
	dP100dt = k1*P000 - (k2a+k3+k4+kret)*P100 + k5a*P1a0 + k6*P101;
	dP0a0dt = k2a*P000 - (k1+k3+k2b+k5a+kret)*P0a0 + k4*P1a0 + k5b*P010 + k6*P0a1;
	dP010dt = k2b*P0a0 -(k1+k3+k5b+kret)*P010 + k4*P110 + k6*P011;
	dP001dt = k3*P000 - (k1+k2a+k6+kret)*P001 + k4*P101 + k5a*P0a1;
	dP1a0dt = k1*P0a0 + k2a*P100 - (k2b+k3+k4+k5a+kret)*P1a0 +k5b*P110+ k6*P1a1;
	dP110dt = k1*P010 + k2b*P1a0 - (k3+k4+k5b+kret+kspli)*P110 + k6*P111;
	dP101dt = k1*P001 + k3*P100 - (k2a+k4+k6+kret+kspli)*P101 + k5a*P1a1;
	dP0a1dt = k2a*P001 + k3*P0a0 - (k1+k2b+k5a+k6+kret)*P0a1 + k4*P1a1 + k5b*P011;
	dP011dt = k2b*P0a1 + k3*P010 - (k1+k5b+k6+kret+kspli)*P011 + k4*P111;
	dP1a1dt = k1*P0a1 + k2a*P101 + k3*P1a0 - (k2b+k4+k5a+k6+kret)*P1a1 + k5b*P111;
	dP111dt = k1*P011 + k2b*P1a1 + k3*P110 - (k4+k5b+k6+kret+(2*kspli))*P111;
	dincldt = kspli*(P111x +Px111) - kincl*incl;
	dskipdt = kspli*P101 - kskip*skip;
	dretdt = kret*(P000+P100+P010+P001+P110+P101+ P011 +P111+P0a0+P1a0+P0a1+P1a1) - (kdr1+kdr2)*ret;
	dP011xdt = kspli*P011 - (k1 + kret)*P011x + k4*P111x;
	dP111xdt = k1*P011x + kspli*P111 - (k4+kret+kspli)*P111x;
	dPx110dt = kspli*P110 - k3*Px110 - kret*Px110 + k6*Px111;
	dPx111dt = k3*Px110 - kret*Px111 + kspli*P111 - (k6+kspli)*Px111;
	dret1dt = kret*(P011x + P111x) - kdr1*ret1;
	dret2dt = kret*(Px110 + Px111) - kdr2*ret2;
	ddegdt = kincl*incl + kskip*skip + kdr1*ret1 + kdr2*ret2 + (kdr1 + kdr2)*ret;
    
	dxdt = [dP000dt;dP100dt;dP0a0dt;dP010dt;dP001dt;dP1a0dt;dP110dt;dP101dt;dP0a1dt;dP011dt;dP1a1dt;dP111dt;dincldt;dskipdt;dretdt;dP011xdt;dP111xdt;dPx110dt;dPx111dt;dret1dt;dret2dt;ddegdt];

end

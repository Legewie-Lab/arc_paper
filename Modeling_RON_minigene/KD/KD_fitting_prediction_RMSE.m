%% KD_fitting_prediction_RMSE.m
% Fit RBP knockdown predictions with the RMSE objective variant.
% Input data and helper functions are expected in the MATLAB working directory.

%% Loading data

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
sol=multik2lsomikronfindercommit_RMSE(alldatwt,indxwt,aa,error_fitpred(1,:));
toc








%% Predict the isoform frequency of mutants using the best-fit parameter

% sol=multik2lsomikronfindercommit(alldatwt,indx,res,num,errorek2);
tic
% val=[];
ksol=[];
for j=1:length(sol)
    val=[];

    for cond=1:length(sol{1})
        y=multiomikroncalculatorcalccommit_RMSE(sol{j,2}(cond),res(indx{j,2},:),alldat{j}(indx{j,1},2:6),0,cond,1);

        val=[val;{y}]; 
%         val(cond,:)={y};
%         cond
    end
    ksol=[ksol;{val}];
end
toc

%%  RMSE prediction score
for j = 1:length(ksol)

    rmse_list = [];

    for cond = 1:15

        % data
        data = alldat{j}(indx{j,1}, 2:6);

        % pred
        pred = ksol{j}{cond};

        %  RMSE
        diff = data - pred;
        rmse = sqrt(mean(diff.^2,'all'));

        rmse_list = [rmse_list; rmse];

    end

    % 2 column
    ksol{j,2} = rmse_list;

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

   solmat(:,i)=ksol{i,2};
end





% %%%%%%%%%%%%%%%%% Figure 5B    %%%%%%%%%%%%%%%%% 
% interleave fitting and prediction columns
[rows, cols] = size(solsd);  
C = zeros(rows, 2 * cols);

for col = 1:cols
    C(:, 2*col - 1) = solsd(:, col);   % fitting RMSE
    C(:, 2*col)     = solmat(:, col);  % prediction RMSE
end

C = log10(C);

% labels
xlabels = {'PRPF6 fitting','PRPF6 prediction',...
           'PUF60 fitting','PUF60 prediction',...
           'HNRNPH1 fitting','HNRNPH1 prediction',...
           'SMU1 fitting','SMU1 prediction',...
           'SRSF2 fitting','SRSF2 prediction',...
           'MCF7 fitting','MCF7 prediction'};

% plot heatmap
figure;
h = heatmap(xlabels, indication, C);

h.Colormap = bone;          
% h.ColorScaling = 'scaled';  
% h.ColorLimits = [0 0.25];   
h.CellLabelFormat = '%.2f';

title('RMSE comparison: Fitting vs Prediction');


% ---  save ---
figure5B_data_RMSE.solsd       = solsd;
figure5B_data_RMSE.solmat      = solmat;
figure5B_data_RMSE.C           = C;
figure5B_data_RMSE.indication  = indication;
figure5B_data_RMSE.xlabels     = xlabels;
figure5B_data_RMSE.metric      = 'RMSE';
figure5B_data_RMSE.transform   = 'log10';
figure5B_data_RMSE.date        = datestr(now);

save('Figure5B_data_RMSE.mat', 'figure5B_data_RMSE');


% ---- plot 1
load("Figure5B_data_RMSE.mat");

figure

ind = [3 2 1 4 5];

for i = 1:5

    col_fit  = 2*ind(i) - 1;
    col_pred = 2*ind(i);

    c = figure5B_data_RMSE.C(:, [col_fit col_pred]);

    subplot(1,5,i)
    h = heatmap(c);

    h.Colormap = bone;
    h.ColorLimits = [-2.2 -1];
    h.ColorScaling = 'scaled';
    h.CellLabelFormat = '%.2f';

end
% --- Plot 2
S = load('Figure5B_data_RMSE.mat');
data = S.figure5B_data_RMSE;

figure;
h = heatmap(data.xlabels, data.indication, data.C);

h.Colormap = bone;
h.CellLabelFormat = '%.2f';

title('RMSE comparison: Fitting vs Prediction');






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


function sol = multik2lsomikronfindercommit_RMSE(wtkd,inde,aa,error)
    
    options=optimset('display','off');
    
    sol = cell(length(wtkd),2);
    
    for j = 1:length(wtkd)
    
        ls = zeros(15,1);
        z  = zeros(15,1);
    
        parfor cond = 1:15
    
            [w,~,residual] = lsqnonlin( ...
                @(b) multiomikroncalculatorcalccommit_RMSE( ...
                    b, aa(inde{j,2},:), ...
                    wtkd{j}(inde{j,1},2:6), ...
                    1, cond, error(inde{j,2},2:6)), ...
                1, 0, 20, options);
    
            rmse = sqrt(mean(residual.^2));
    
            ls(cond) = rmse;
            z(cond)  = w;
    
        end
    
        sol{j,1} = ls;
        sol{j,2} = z;
    
    end
end



function [y,o] = multiomikroncalculatorcalccommit_RMSE(b,aa,output,fit,cond,error)
    
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
        y = (Y1-output);
    
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
